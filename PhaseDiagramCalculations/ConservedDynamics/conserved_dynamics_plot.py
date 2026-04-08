from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def compute_surface_coordinates(mesh, phi, radial_scale=-1.0):
    x3d, y3d, z3d = mesh.faceCenters
    x3d = np.asarray(x3d)
    y3d = np.asarray(y3d)
    z3d = np.asarray(z3d)
    values = np.asarray(phi.arithmeticFaceValue)

    azimuth = np.arctan2(y3d, x3d)
    radius = np.sqrt(x3d ** 2 + y3d ** 2 + z3d ** 2)
    polar = np.arccos(z3d / radius)

    displaced_radius = radius + radial_scale * values
    x_coords = displaced_radius * np.cos(azimuth) * np.sin(polar)
    y_coords = displaced_radius * np.sin(azimuth) * np.sin(polar)
    z_coords = displaced_radius * np.cos(polar)

    value_min = values.min()
    value_max = values.max()
    if np.isclose(value_max, value_min):
        normalized = np.zeros_like(values)
    else:
        normalized = (values - value_min) / (value_max - value_min)

    return x_coords, y_coords, z_coords, normalized


def render_surface_plot(
    x_coords,
    y_coords,
    z_coords,
    values,
    output_path=None,
    show=False,
    title=None,
    point_size=8.0,
    cmap="YlOrBr",
):
    import matplotlib.pyplot as plt

    figure = plt.figure(figsize=(8, 8))
    axis = figure.add_subplot(111, projection="3d")
    scatter = axis.scatter(
        x_coords,
        y_coords,
        z_coords,
        c=values,
        cmap=cmap,
        s=point_size,
        linewidths=0,
    )
    axis.set_axis_off()
    axis.set_box_aspect((1.0, 1.0, 1.0))

    if title:
        axis.set_title(title)

    figure.colorbar(scatter, ax=axis, shrink=0.65, pad=0.02, label="Normalized field")

    if output_path:
        figure.savefig(output_path, dpi=300, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(figure)

    return output_path


def render_solution(
    mesh,
    phi,
    output_path=None,
    show=False,
    title=None,
    radial_scale=-1.0,
    point_size=8.0,
):
    x_coords, y_coords, z_coords, values = compute_surface_coordinates(mesh, phi, radial_scale=radial_scale)
    return render_surface_plot(
        x_coords,
        y_coords,
        z_coords,
        values,
        output_path=output_path,
        show=show,
        title=title,
        point_size=point_size,
    )


def render_dump(
    input_path,
    output_path=None,
    show=False,
    title=None,
    radial_scale=-1.0,
    point_size=8.0,
):
    from fipy.tools import dump

    mesh, phi = dump.read(str(input_path))
    return render_solution(
        mesh,
        phi,
        output_path=output_path,
        show=show,
        title=title,
        radial_scale=radial_scale,
        point_size=point_size,
    )


def default_output_path(input_path):
    input_path = Path(input_path)
    return input_path.with_name(input_path.name + ".png")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Render a conserved-dynamics FiPy dump without mayavi."
    )
    parser.add_argument("input_path", help="Path to a FiPy dump file such as q0d5, q1d5, or q2.")
    parser.add_argument(
        "--output",
        help="Output image path. Defaults to <input>.png when --show is not used.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the figure interactively after rendering.",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional plot title.",
    )
    parser.add_argument(
        "--radial-scale",
        type=float,
        default=-1.0,
        help="Radial displacement factor applied to the field values.",
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=8.0,
        help="Scatter point size used for the rendered surface.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    output_path = args.output

    if output_path is None and not args.show:
        output_path = default_output_path(args.input_path)

    render_dump(
        args.input_path,
        output_path=output_path,
        show=args.show,
        title=args.title,
        radial_scale=args.radial_scale,
        point_size=args.point_size,
    )

    if output_path:
        print(output_path)


if __name__ == "__main__":
    main()
