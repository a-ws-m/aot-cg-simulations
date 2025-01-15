"""Tools for analysing bilayer thickness and area per lipid."""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

try:
    from functools import cached_property
except ImportError:
    from cached_property import cached_property

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd
import pytim
import pyvista as pv
import seaborn as sns
from MDAnalysis.analysis.distances import distance_array
from pytim.datafiles import CHARMM27_TOP, pytim_data
from scipy.fft import rfft, rfftfreq
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
from tqdm import trange


def self_distances(atomgroup: mda.AtomGroup, max_dist: float = 6.0):
    """Calculate the self radial distribution function for a given atomgroup."""
    dist_arr = distance_array(atomgroup, atomgroup, box=atomgroup.dimensions)
    one_pair_each = np.triu(dist_arr, k=1).flatten()
    one_pair_each = one_pair_each[one_pair_each > 0.0]
    one_pair_each = one_pair_each[one_pair_each <= max_dist]
    return one_pair_each


def partition_xy(atomgroup: mda.AtomGroup, n_bins: int = 20):
    """Partition the xy plane into (n_bins*n_bins) bin indexes for each atom."""
    bin_edges = np.arange(0, atomgroup.dimensions[0], atomgroup.dimensions[0] / n_bins)
    return np.digitize(atomgroup.positions[:, :2], bin_edges)

def rt_distance(surf: pv.PolyData, return_surf: bool = False) -> "float":
    """Calculate the ray traced distance between the upper and lower surface."""
    split = surf.split_bodies()
    split = split.as_polydata_blocks()

    surf0_normals = split[0].compute_normals(point_normals=True, cell_normals=False, auto_orient_normals=True)
    surf0_normals["distances"] = np.empty(surf0_normals.n_points)

    for i in range(surf0_normals.n_points):
        p = surf0_normals.points[i]
        vec = surf0_normals["Normals"][i] * surf0_normals.length
        p0 = p - vec
        p1 = p + vec
        ip, ic = split[1].ray_trace(p0, p1, first_point=True)
        dist = np.sqrt(np.sum((ip - p) ** 2))
        surf0_normals["distances"][i] = dist
    
    # Replace zeros with nans
    mask = surf0_normals["distances"] == 0
    surf0_normals["distances"][mask] = np.nan

    if return_surf:
        return surf0_normals, split[1]
    else:
        return np.nanmean(surf0_normals["distances"])

def nn_distance(surf: pv.PolyData, return_surf: bool = False) -> "float":
    """Calculate the nearest-neighbour distance between the upper and lower surface."""
    split = surf.split_bodies()
    split = split.as_polydata_blocks()

    if len(split) != 2 or np.abs(split[0].n_points - split[1].n_points) > 1000:
        raise ValueError("Surfaces are not separated properly.")

    tree = KDTree(split[1].points)
    d_kdtree, idx = tree.query(split[0].points)
    split[0]["distances"] = d_kdtree

    if return_surf:
        return split[0], split[1]
    
    else:
        return np.mean(split[0]["distances"])

def sep_willard_chandler(
    group: mda.AtomGroup,
    radii_dict: dict = pytim_data.vdwradii(CHARMM27_TOP),
) -> "float":
    """Plot the separated Willard-Chandler surface for a given group of atoms."""
    u = group.universe

    # Get the surface
    wc = pytim.WillardChandler(
        u,
        group=group,
        alpha=3.0,
        mesh=2.0,
        radii_dict=radii_dict,
        autoassign=False,
        density_cutoff=0.016,
        # centered=True,
    )

    # radius, _, _, _ = pytim.utilities.fit_sphere(wc.triangulated_surface[0])

    # converting PyTim to PyVista surface
    verts = wc.triangulated_surface[0]
    faces = wc.triangulated_surface[1]
    threes = 3 * np.ones((faces.shape[0], 1), dtype=int)
    faces = np.concatenate((threes, faces), axis=1)
    poly = pv.PolyData(verts, faces)
    
    # surf0, surf1 = nn_distance(poly, return_surf=True)

    # p = pv.Plotter()
    # p.add_mesh(surf0, scalars="distances", smooth_shading=True)
    # p.add_mesh(surf1, color=True, opacity=0.75, smooth_shading=True)
    # p.show()

def willard_chandler(
    group: mda.AtomGroup,
    radii_dict: dict = pytim_data.vdwradii(CHARMM27_TOP),
    save: Optional[str] = None,
    save_thickness: Optional[float] = None,
) -> "tuple[float, float]":
    """Calculate the Willard-Chandler surface for a given group of atoms and associated properties."""
    u = group.universe

    # Get the surface
    wc = pytim.WillardChandler(
        u,
        group=group,
        alpha=4.0,
        mesh=2.0,
        radii_dict=radii_dict,
        autoassign=False,
        density_cutoff_ratio=0.33,
        # centered=True,
    )

    if save is not None:
        wc.writecube(save + ".cube", group=group, normalize=False)

    # radius, _, _, _ = pytim.utilities.fit_sphere(wc.triangulated_surface[0])

    # converting PyTim to PyVista surface
    verts = wc.triangulated_surface[0]
    faces = wc.triangulated_surface[1]
    threes = 3 * np.ones((faces.shape[0], 1), dtype=int)
    faces = np.concatenate((threes, faces), axis=1)
    poly = pv.PolyData(verts, faces)

    # Get actual surface volume
    area = poly.area

    # Get thickness
    thickness = nn_distance(poly)

    if save_thickness is not None and thickness >= save_thickness:
        surf0, surf1 = nn_distance(poly, return_surf=True)
        wc.writepdb(f"{thickness:.2f}-nm-thickness.pdb")
        surf0.save(f"{thickness:.2f}-nm-thickness-upper.stl", binary=False)
        surf1.save(f"{thickness:.2f}-nm-thickness-lower.stl", binary=False)
        print(f"Thickness = {thickness:.2f} A, lower points = {len(surf1.points)}, upper points = {len(surf0.points)}")

    return area, thickness


def binned_thickness(atomgroup: mda.AtomGroup, n_bins: int = 20):
    """Calculate the thickness distribution for each bin."""
    indexed = partition_xy(atomgroup, n_bins)
    thicknesses = []
    for i in range(n_bins):
        for j in range(n_bins):
            this_index = np.where(indexed == [i + 1, j + 1])
            if not np.any(this_index):
                raise ValueError("No atoms found in this bin.")

            these_z_vals = atomgroup.positions[this_index, 2]
            avg_z_val = np.mean(these_z_vals)
            upper_layer = these_z_vals[these_z_vals >= avg_z_val]
            lower_layer = these_z_vals[these_z_vals < avg_z_val]

            # thicknesses.extend(cdist(upper_layer[:, np.newaxis], lower_layer[:, np.newaxis], "cityblock").flatten())
            thicknesses.append(np.mean(upper_layer) - np.mean(lower_layer))

    return thicknesses


def binned_mean_z(atomgroup: mda.AtomGroup, n_bins: int = 20):
    """Compute the average z values for a given atomgroup."""
    indexed = partition_xy(atomgroup, n_bins)
    avg_z_vals = np.zeros((n_bins, n_bins))
    for i in range(n_bins):
        for j in range(n_bins):
            this_index = np.where(indexed == [i + 1, j + 1])
            if not np.any(this_index):
                raise ValueError("No atoms found in this bin.")
            avg_z_vals[i, j] = np.mean(atomgroup.positions[this_index, 2])

    return avg_z_vals


class Bilayer:
    """Class for analysing bilayer thickness and area per lipid."""

    def __init__(
        self,
        gro_file,
        xtc_file,
        name,
        head_selection="name SU",
        tail_selection="name C[LR]*",
    ):
        self.gro_file = gro_file
        self.xtc_file = xtc_file
        self.name = name
        self.selection = head_selection
        self.u = mda.Universe(gro_file, xtc_file)
        self.head_atoms = self.u.select_atoms(head_selection)
        self.tail_atoms = self.u.select_atoms(tail_selection)
        if not isinstance(self.head_atoms, mda.AtomGroup):
            raise ValueError("No atoms were selected. Check your selection string.")

        self.aot: mda.AtomGroup = self.u.select_atoms("resname AOT")
        self.water = self.u.select_atoms("resname W")

        aot_wt = len(self.aot) * 444.6
        water_wt = len(self.water) * 18.01528
        self.wt_percent = 100 * (aot_wt) / (aot_wt + water_wt)

    @cached_property
    def vdwradii(self) -> "dict[str, float]":
        """Determine the van der Waals radii for the atoms in the system."""
        radii = dict()
        for atom in self.u.atoms:
            # See Section 8.2 of
            # https://cgmartini.nl/docs/tutorials/Martini3/Small_Molecule_Parametrization/#molecular-volume-and-shape
            atom_type = atom.type[0]
            if atom_type == "S":
                radius = 0.230
            elif atom_type == "T":
                radius = 0.191
            else:
                radius = 0.264

            radii[atom.name] = radius

        return radii

    def _get_rdfs(self, selection: mda.AtomGroup, n_steps=20) -> pd.DataFrame:
        """Get the RDFs for a given selection."""
        step_size = self.u.trajectory.n_frames // n_steps
        plot_data = {"Time (ns)": [], "RDF": []}
        for i in trange(
            0, self.u.trajectory.n_frames, step_size, desc="Calculating RDFs"
        ):
            self.u.trajectory[i]
            self_dists = self_distances(selection)
            plot_data["RDF"].extend(list(self_dists))
            plot_data["Time (ns)"].extend(
                [self.u.trajectory.time / 1000] * len(self_dists)
            )

        return pd.DataFrame(plot_data)

    def get_head_rdfs(self, n_steps=10) -> pd.DataFrame:
        """Get the RDFs for the headgroup atoms."""
        return self._get_rdfs(self.head_atoms, n_steps)

    def get_tail_rdfs(self, n_steps=10) -> pd.DataFrame:
        """Get the RDFs for the headgroup atoms."""
        return self._get_rdfs(self.tail_atoms, n_steps)

    def new_thickness_dist(self, start_frame=0, end_frame=None, n_bins=20):
        """Calculate the thickness distribution for the bilayer."""
        if end_frame is None:
            end_frame = self.u.trajectory.n_frames

        data = {r"Thickness ($\AA$)": [], "Time (ns)": []}
        for i in trange(start_frame, end_frame, desc="Calculating thicknesses"):
            self.u.trajectory[i]
            current_thicknesses = binned_thickness(self.head_atoms, n_bins)
            data[r"Thickness ($\AA$)"].extend(current_thicknesses)
            data["Time (ns)"].extend(
                [self.u.trajectory.time / 1000] * len(current_thicknesses)
            )

        return data

    # def curved_area_per_headgroup(
    #     self, start_frame=0, end_frame=None, step=10, save_step: Optional[int] = None
    # ):
    #     """Calculate the area per headgroup for a curved bilayer."""
    #     if end_frame is None:
    #         end_frame = self.u.trajectory.n_frames

    #     plot_data = {"Time (ns)": [], r"Area per headgroup ($\AA^2$)": []}
    #     for i in trange(
    #         start_frame,
    #         end_frame,
    #         step,
    #         desc="Calculating area per headgroup for curved bilayer",
    #     ):
    #         self.u.trajectory[i]

    #         save = f"{self.name}-interface-{i}.cube" if i == save_step else None
    #         area = willard_chandler(self.aot, radii_dict=self.vdwradii, save=save)

    #         plot_data["Time (ns)"].append(self.u.trajectory.time / 1000)
    #         plot_data[r"Area per headgroup ($\AA^2$)"].append(
    #             area / self.aot.n_residues
    #         )

    #     return plot_data

    # def plot_curved_area_per_headgroup(
    #     self, output_file, start_frame=0, end_frame=None, step=10
    # ):
    #     """Plot the area per headgroup for a curved bilayer."""
    #     plot_data = self.curved_area_per_headgroup(start_frame, end_frame, step)
    #     plt.figure()
    #     sns.scatterplot(
    #         data=plot_data,
    #         x="Time (ns)",
    #         y=r"Area per headgroup ($\AA^2$)",
    #     )
    #     plt.axvline(64, color="k", linestyle="--", zorder=1.5, alpha=0.5)
    #     plt.tight_layout()
    #     plt.savefig(output_file)

    def wc_properties(
        self, start_frame=0, end_frame=None, step=10, save_step: Optional[int] = None, save_thickness: Optional[float] = None
    ):
        """Calculate the area per headgroup and thickness for a curved bilayer."""
        if end_frame is None:
            end_frame = self.u.trajectory.n_frames

        plot_data = {"Time (ns)": [], r"Area per headgroup ($\AA^2$)": [], r"Thickness ($\AA$)": []}
        for i in trange(
            start_frame,
            end_frame,
            step,
            desc="Calculating Willard-Chandler properties for curved bilayer",
        ):
            self.u.trajectory[i]

            save = f"{self.name}-interface-{i}" if i == save_step else None
            try:
                area, thickness = willard_chandler(self.aot, radii_dict=self.vdwradii, save=save, save_thickness=save_thickness)
            except ValueError:
                # Can't separate the surfaces
                print(f"Couldn't separate surfaces at {self.u.trajectory.time} ps")
                continue

            plot_data["Time (ns)"].append(self.u.trajectory.time / 1000)
            plot_data[r"Area per headgroup ($\AA^2$)"].append(
                area / self.aot.n_residues
            )
            plot_data[r"Thickness ($\AA$)"].append(thickness)

        return plot_data

    def plot_new_thickness_trend(
        self, output_file, start_frame=0, end_frame=None, n_bins=20
    ):
        """Plot the bilayer thicknesses over time."""
        plot_data = self.new_thickness_dist(start_frame, end_frame, n_bins)
        plt.figure()
        sns.lineplot(
            data=plot_data,
            x="Time (ns)",
            y=r"Thickness ($\AA$)",
            errorbar="sd",
        )
        plt.axvline(20, color="k", linestyle="--", zorder=1.5, alpha=0.5)
        plt.tight_layout()
        plt.savefig(output_file)

    def plot_new_thickness_dist(
        self, output_file, start_frame=0, end_frame=None, n_bins=20
    ):
        """Plot the distribution of bilayer thicknesses."""
        plot_data = self.new_thickness_dist(start_frame, end_frame, n_bins)
        plt.figure()
        sns.violinplot(
            plot_data,
            x=r"Thickness ($\AA$)",
        )
        plt.tight_layout()
        plt.savefig(output_file)

    def plot_oop(self, output_file, n_steps=10):
        step_size = self.u.trajectory.n_frames // n_steps
        plot_data = {"Time (ns)": [], r"Displacement ($\AA$)": []}
        for i in trange(
            0, self.u.trajectory.n_frames, step_size, desc="Calculating OOP density"
        ):
            self.u.trajectory[i]
            fitted = PCA(n_components=3).fit_transform(self.head_atoms.positions)
            oop_coords = fitted[:, 2]
            oop_coords -= np.mean(oop_coords)
            plot_data[r"Displacement ($\AA$)"].extend(list(oop_coords))
            plot_data["Time (ns)"].extend(
                [self.u.trajectory.time / 1000] * len(self.head_atoms)
            )

        plt.figure()
        sns.violinplot(
            plot_data,
            x="Time (ns)",
            y=r"Displacement ($\AA$)",
            inner=None,
            density_norm="area",
            common_norm=True,
            bw_adjust=0.2,
            native_scale=True,
        )
        plt.tight_layout()
        plt.savefig(output_file, transparent=True)

    def get_xy_area(self):
        """Calculate the xy area of the bilayer."""
        data = {"Time (ns)": [], r"$xy$-Area ($\AA^2$)": []}
        for i in trange(self.u.trajectory.n_frames, desc="Calculating xy area"):
            self.u.trajectory[i]
            data[r"$xy$-Area ($\AA^2$)"].append(np.prod(self.u.dimensions[:2]))
            data["Time (ns)"].append(self.u.trajectory.time / 1000)

        return data

    def get_thickness_dist(self, n_steps=10):
        step_size = self.u.trajectory.n_frames // n_steps

        plot_data = {"Time (ns)": [], r"Thickness ($\AA$)": []}

        for i in trange(
            0, self.u.trajectory.n_frames, step_size, desc="Calculating thicknesses"
        ):
            self.u.trajectory[i]
            fitted = PCA(n_components=3).fit_transform(self.head_atoms.positions)
            oop_coords = fitted[:, 2]
            oop_coords -= np.mean(oop_coords)

            top_layer = oop_coords > 0
            thicknesses = cdist(
                oop_coords[top_layer, np.newaxis],
                oop_coords[~top_layer, np.newaxis],
                "cityblock",
            )
            plot_data[r"Thickness ($\AA$)"].extend(list(thicknesses.flatten()))
            plot_data["Time (ns)"].extend(
                [self.u.trajectory.time / 1000] * len(thicknesses.flatten())
            )

        return plot_data

    def plot_thickness_dist(self, output_file, n_steps=10, **savefig_kwargs):
        """Plot the distribution of bilayer thicknesses."""
        plot_data = self.get_thickness_dist(n_steps)
        plt.figure()
        sns.violinplot(
            plot_data,
            x="Time (ns)",
            y=r"Thickness ($\AA$)",
            density_norm="area",
            inner="quart",
            common_norm=True,
            bw_adjust=0.2,
            native_scale=True,
        )
        plt.tight_layout()
        plt.savefig(output_file, **savefig_kwargs)

    def get_area_per_headgroup(self):
        plot_data = {"Time (ns)": [], r"Area per headgroup ($\AA^2$)": []}
        for i in trange(
            self.u.trajectory.n_frames, desc="Calculating area per headgroup"
        ):
            self.u.trajectory[i]
            ip_area = np.prod(self.u.dimensions[:2])
            area_per_hg = 2 * ip_area / len(self.head_atoms)

            plot_data["Time (ns)"].append(self.u.trajectory.time / 1000)
            plot_data[r"Area per headgroup ($\AA^2$)"].append(area_per_hg)

        return plot_data

    def plot_area_per_headgroup(self, output_file, **savefig_kwargs):
        plot_data = self.get_area_per_headgroup()

        plt.figure()
        g = sns.relplot(
            data=plot_data,
            x="Time (ns)",
            y=r"Area per headgroup ($\AA^2$)",
            kind="line",
        )
        g.tight_layout()
        g.savefig(output_file, **savefig_kwargs)

    def plot_head_rdfs(self, output_file, n_steps=10, **savefig_kwargs):
        """Plot the RDFs for the headgroup atoms."""
        rdf_data = self.get_head_rdfs(n_steps)
        plt.figure()
        sns.violinplot(
            rdf_data,
            x="Time (ns)",
            y="RDF",
            inner=None,
            density_norm="area",
            common_norm=True,
            bw_adjust=0.2,
            native_scale=True,
        )
        plt.savefig(output_file, **savefig_kwargs)

    def plot_tail_rdfs(self, output_file, n_steps=10, **savefig_kwargs):
        """Plot the RDFs for the headgroup atoms."""
        rdf_data = self.get_head_rdfs(n_steps)
        plt.figure()
        sns.violinplot(
            rdf_data,
            x="Time (ns)",
            y="RDF",
            inner=None,
            density_norm="area",
            common_norm=True,
            bw_adjust=0.2,
            native_scale=True,
        )
        plt.savefig(output_file, **savefig_kwargs)

    def plot_z_fft(
        self, output_file, n_bins: int = 20, n_steps=10, spacing=5, **savefig_kwargs
    ):
        """Calculate the 2D fft"""
        fig = plt.figure(figsize=(8, 8), facecolor="black")
        ax = plt.subplot(frameon=False)

        step_size = self.u.trajectory.n_frames // n_steps

        frame_ffts = []
        frame_freqs = []
        all_mean_zs = []
        for i in trange(
            0, self.u.trajectory.n_frames, step_size, desc="Calculating ffts"
        ):
            self.u.trajectory[i]
            mean_z_vals = binned_mean_z(self.u.select_atoms("resname AOT"), n_bins)
            all_mean_zs.append(mean_z_vals)
            this_fft = []
            this_freq = []
            for j in range(n_bins):
                this_fft.append(np.abs(rfft(mean_z_vals[j])))
                this_freq.append(rfftfreq(n_bins, d=self.u.dimensions[0] / n_bins))
            frame_ffts.append(np.array(this_fft))
            frame_freqs.append(np.array(this_freq))

        # plt.plot(frame_freqs[0][0], frame_ffts[0][0])
        # plt.xlabel("Frequency")
        # plt.savefig("example-f.png")

        # plt.figure()
        # plt.plot(all_mean_zs[0][0])
        # plt.xlabel("Distance")
        # plt.savefig("example-x.png")

        lines = []
        for i in range(len(frame_ffts[0])):
            (line,) = ax.plot(
                frame_freqs[0][i], np.log(frame_ffts[0][i]) + i * spacing, color="w"
            )
            lines.append(line)

        ax.set_yticks([])
        ax.tick_params(axis="x", colors="white")
        ax.set_xlabel("Frequency (1/Angstrom)")
        ax.xaxis.label.set_color("white")
        ax.set_xlim(right=0.04)

        def update(frame):
            for idx, line in enumerate(lines):
                line.set_xdata(frame_freqs[frame][idx])
                line.set_ydata(np.log(frame_ffts[frame][idx]) + idx * spacing)

            ax.tick_params(axis="x", colors="white")
            ax.set_xlabel("Frequency (1/Angstrom)")
            ax.xaxis.label.set_color("white")

            return lines

        ani = animation.FuncAnimation(fig, update, frames=n_steps, blit=True)
        ani.save(
            output_file, writer="ffmpeg", fps=2, dpi=300, savefig_kwargs=savefig_kwargs
        )


@dataclass
class ModelFiles:
    """Contains locations of files for specific model."""

    name: str
    friendly_name: str
    xtc_name: str
    gro_name: str
    head_name: str
    tail_name: str

    def __post_init__(self):
        self.dir = Path(__file__).parent / f"{self.name}-bilayer-new"
        self.xtc = self.dir / self.xtc_name
        self.gro = self.dir / self.gro_name

        self.head_selection = f"name {self.head_name}"
        self.tail_selection = f"name {self.tail_name}"

        for file in [self.dir, self.xtc, self.gro]:
            if not file.exists():
                raise FileNotFoundError(f"{file} not found.")


class BilayerPlotter:
    """Compare plots for several bilayer systems."""

    def __init__(self, models: list[ModelFiles]):
        self.models = models
        self.bilayers = {
            m.friendly_name: Bilayer(
                m.gro, m.xtc, m.friendly_name, m.head_selection, m.tail_selection
            )
            for m in models
        }

    def plot_xy_area(self, output_file, transparent=False):
        df = pd.DataFrame()
        for name, bilayer in self.bilayers.items():
            areas = bilayer.get_xy_area()
            areas["Model"] = name
            df = pd.concat([df, pd.DataFrame(areas)], ignore_index=True)

            these_areas = areas[r"$xy$-Area ($\AA^2$)"]
            print(f"{name} mean xy area = {np.mean(these_areas):.2f}")

        plt.figure()
        g = sns.relplot(
            data=df,
            x="Time (ns)",
            y=r"$xy$-Area ($\AA^2$)",
            hue="Model",
            kind="line",
        )
        g.tight_layout()
        g.savefig(output_file, transparent=transparent)

    def plot_new_thickness(
        self,
        output_file,
        n_bins=10,
        start_frame: int = 0,
        end_frame=None,
        transparent=False,
    ):
        df = pd.DataFrame()
        for name, bilayer in self.bilayers.items():
            thicknesses = bilayer.new_thickness_dist(start_frame, end_frame, n_bins)
            thicknesses["Model"] = name
            df = pd.concat([df, pd.DataFrame(thicknesses)], ignore_index=True)

            thickness_arr = thicknesses[r"Thickness ($\AA$)"]
            print(f"{name} mean thickness = {np.mean(thickness_arr):.2f}")

        plt.figure()
        g = sns.catplot(
            df,
            x=r"Thickness ($\AA$)",
            y="Model",
            hue="Model",
            kind="violin",
            density_norm="width",
            fill=False,
            common_norm=True,
        )
        plt.axvline(20, color="k", linestyle="--", zorder=1.5, alpha=0.5)
        g.tight_layout()
        g.savefig(output_file, transparent=transparent)

    # def plot_curved_area_per_headgroup(
    #     self,
    #     output_file,
    #     start_frame=0,
    #     end_frame=None,
    #     step=10,
    #     transparent: bool = False,
    #     save_step: Optional[int] = None,
    # ):
    #     df = pd.DataFrame()
    #     for name, bilayer in self.bilayers.items():
    #         areas = bilayer.curved_area_per_headgroup(
    #             start_frame, end_frame, step, save_step=save_step
    #         )
    #         areas["Model"] = name
    #         df = pd.concat([df, pd.DataFrame(areas)], ignore_index=True)

    #         aphg = areas[r"Area per headgroup ($\AA^2$)"]
    #         print(f"{name} mean area per headgroup = {np.mean(aphg):.2f}")

    #     plt.figure()
    #     g = sns.catplot(
    #         df,
    #         x=r"Area per headgroup ($\AA^2$)",
    #         y="Model",
    #         hue="Model",
    #         kind="violin",
    #         density_norm="width",
    #         fill=False,
    #         common_norm=True,
    #     )
    #     plt.axvline(64, color="k", linestyle="--", zorder=1.5, alpha=0.5)
    #     g.tight_layout()
    #     g.savefig(output_file, transparent=transparent)

    def plot_wc_properties(
        self,
        area_file,
        thickness_file,
        start_frame=0,
        end_frame=None,
        step=10,
        transparent: bool = False,
        save_step: Optional[int] = None,
    ):
        df = pd.DataFrame()
        for name, bilayer in self.bilayers.items():
            wc_props = bilayer.wc_properties(
                start_frame, end_frame, step, save_step=save_step
            )
            wc_props["Model"] = name
            df = pd.concat([df, pd.DataFrame(wc_props)], ignore_index=True)

            aphg = wc_props[r"Area per headgroup ($\AA^2$)"]
            print(f"{name} mean area per headgroup = {np.mean(aphg):.2f}")
            thickness = wc_props[r"Thickness ($\AA$)"]
            print(f"{name} mean thickness = {np.mean(thickness):.2f}")

        plt.figure()
        g = sns.catplot(
            df,
            x=r"Area per headgroup ($\AA^2$)",
            y="Model",
            hue="Model",
            kind="violin",
            density_norm="width",
            fill=False,
            common_norm=True,
        )
        plt.axvline(64, color="k", linestyle="--", zorder=1.5, alpha=0.5)
        g.tight_layout()
        g.savefig(area_file, transparent=transparent)

        plt.figure()
        g = sns.catplot(
            df,
            x=r"Thickness ($\AA$)",
            y="Model",
            hue="Model",
            kind="violin",
            density_norm="width",
            fill=False,
            common_norm=True,
        )
        plt.axvline(20, color="k", linestyle="--", zorder=1.5, alpha=0.5)
        g.tight_layout()
        g.savefig(thickness_file, transparent=transparent)

    def plot_thickness(self, output_file, n_steps=10, transparent=False):
        df = pd.DataFrame()
        for name, bilayer in self.bilayers.items():
            thicknesses = bilayer.get_thickness_dist(n_steps)
            thicknesses["Model"] = name
            df = pd.concat([df, pd.DataFrame(thicknesses)], ignore_index=True)

        plt.figure()
        g = sns.catplot(
            data=df,
            x="Time (ns)",
            y=r"Thickness ($\AA$)",
            hue="Model",
            kind="violin",
            inner="quart",
            # split=True,
            density_norm="width",
            common_norm=True,
            bw_adjust=0.2,
            native_scale=True,
            aspect=1.5,
        )
        g.tight_layout()
        g.savefig(output_file, transparent=transparent)

    def plot_areas_per_hg(
        self, output_file, target: Optional[float] = None, transparent=False
    ):
        df = pd.DataFrame()
        for name, bilayer in self.bilayers.items():
            areas = bilayer.get_area_per_headgroup()
            areas["Model"] = name
            areas_df = pd.DataFrame(areas)
            areas_df["Expanding average"] = (
                areas_df[r"Area per headgroup ($\AA^2$)"].expanding().mean()
            )
            df = pd.concat([df, areas_df], ignore_index=True)

        df[r"Area / headgroup ($\AA^2$)"] = df[r"Area per headgroup ($\AA^2$)"]
        plt.figure()
        g = sns.relplot(
            data=df,
            x="Time (ns)",
            y="Expanding average",
            hue="Model",
            kind="line",
            linestyle="-",
        )
        g.map_dataframe(
            sns.scatterplot,
            data=df,
            x="Time (ns)",
            y=r"Area / headgroup ($\AA^2$)",
            hue="Model",
            marker=".",
        )
        if target:
            plt.axhline(target, color="k", linestyle=":")
        # g.tight_layout()
        g.savefig(output_file, transparent=transparent)


if __name__ == "__main__":
    sns.set_theme(context="paper", style="white")
    zeta = ModelFiles("zeta", "Finest", "3-run.xtc", "no-water.pdb", "SULF", "T*")
    coarse = ModelFiles(
        "coarse-alpha",
        "Coarsest",
        "3-run.xtc",
        "no-water.pdb",
        "SU",
        "C[MB][LR]",
    )
    mixed = ModelFiles(
        "mixed-alpha", "Mixed", "3-run.xtc", "no-water.pdb", "SU", "C[LR]*"
    )

    fine_bilayer = Bilayer(zeta.gro, zeta.xtc, zeta.friendly_name, zeta.head_selection, zeta.tail_selection)
    # coarse_bilayer = Bilayer(coarse.gro, coarse.xtc, coarse.friendly_name, coarse.head_selection, coarse.tail_selection)
    # coarse_bilayer.wc_properties(150, save_thickness=30)
    # fine_bilayer.u.trajectory[150]
    # sep_willard_chandler(fine_bilayer.aot, fine_bilayer.vdwradii)
    # print(f"{fine_bilayer.u.trajectory.ts.dt=}")
    # print(f"{3e+4 / fine_bilayer.u.trajectory.ts.dt=:.2f}")

    model_files = [zeta, mixed, coarse]
    bilayer_plotter = BilayerPlotter(model_files)
    bilayer_plotter.plot_wc_properties(
        "wc-area.pdf",
        "wc-thickness.pdf",
        start_frame=150,
        step=5,
        transparent=True,
        save_step=int(150e3 / fine_bilayer.u.trajectory.ts.dt),
    )
    # bilayer_plotter.plot_new_thickness(
    #     "new-thickness-comparison.pdf", start_frame=150, transparent=True
    # )
    # bilayer_plotter.plot_curved_area_per_headgroup(
    #     "curved-area-comparison.pdf",
    #     start_frame=150,
    #     step=5,
    #     transparent=True,
    #     save_step=int(150e3 / fine_bilayer.u.trajectory.ts.dt),
    # )
    # bilayer_plotter.plot_xy_area("xy-area-comparison.pdf", transparent=True)