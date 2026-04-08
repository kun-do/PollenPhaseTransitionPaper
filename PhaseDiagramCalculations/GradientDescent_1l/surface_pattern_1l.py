import argparse

import numpy as np

try:
    from scipy.special import sph_harm as _spherical_harmonic_impl

    def spherical_harmonic(m_value, ell, theta, phi):
        return _spherical_harmonic_impl(m_value, ell, theta, phi)

except ImportError:
    from scipy.special import sph_harm_y as _spherical_harmonic_impl

    def spherical_harmonic(m_value, ell, theta, phi):
        return _spherical_harmonic_impl(ell, m_value, phi, theta)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot a one-l equilibrium surface pattern on a sphere."
    )
    parser.add_argument("cms", help="Path to endcms.txt or another coefficient file.")
    parser.add_argument("parameters", help="Path to parameters.txt.")
    parser.add_argument("--output", help="Optional output image path.")
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the plot interactively even when --output is provided.",
    )
    return parser.parse_args()


def load_coefficients(path):
    real_parts = []
    imaginary_parts = []

    with open(path, "r") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 2:
                continue
            real_parts.append(float(fields[0]))
            imaginary_parts.append(float(fields[1]))

    return real_parts, imaginary_parts


def load_parameters(path):
    parameters = []

    with open(path, "r") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 2:
                continue
            parameters.append(float(fields[1]))
            if len(parameters) >= 4:
                break

    with open(path, "r") as handle:
        lines = handle.readlines()

    h_value = lines[-1].split()[1] if lines else "unknown"
    parameters.append(h_value)
    return parameters


def build_color_field(real_parts, imaginary_parts, theta, phi):
    ell = len(real_parts) - 1
    colors = (real_parts[0] * spherical_harmonic(0, ell, theta, phi)).real

    for m_value in range(1, ell + 1):
        coefficient = complex(real_parts[m_value], imaginary_parts[m_value])
        colors += (coefficient * spherical_harmonic(m_value, ell, theta, phi)).real
        colors += (
            coefficient.conjugate() * pow(-1, m_value) * spherical_harmonic(-m_value, ell, theta, phi)
        ).real

    color_min = colors.min()
    color_max = colors.max()
    if np.isclose(color_max, color_min):
        return np.zeros_like(colors)
    return (colors - color_min) / (color_max - color_min)


def plot_pattern(real_parts, imaginary_parts, parameters, output_path=None, show=False):
    import matplotlib as mpl

    if output_path and not show:
        mpl.use("Agg")

    import matplotlib.pyplot as plt
    from matplotlib import cm
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    phi = np.linspace(0, np.pi, 150)
    theta = np.linspace(0, 2 * np.pi, 150)
    phi, theta = np.meshgrid(phi, theta)

    x_coords = np.sin(phi) * np.cos(theta)
    y_coords = np.sin(phi) * np.sin(theta)
    z_coords = np.cos(phi)

    colors = build_color_field(real_parts, imaginary_parts, theta, phi)

    figure = plt.figure(figsize=(10, 5))
    axis = figure.add_subplot(1, 2, 1, projection="3d")
    axis.plot_surface(x_coords, y_coords, z_coords, rstride=1, cstride=1, facecolors=cm.seismic(colors))
    axis.view_init(45, 0)
    axis.set_axis_off()

    axis = figure.add_subplot(1, 2, 2, projection="3d")
    axis.plot_surface(x_coords, y_coords, z_coords, rstride=1, cstride=1, facecolors=cm.seismic(colors))
    axis.view_init(45, 180)
    axis.set_axis_off()

    plt.title(
        "Parameters are: el_not=%s, tau=%s, lambda3=%s, lambda4=%s\nMinimum Hamiltonian value: %s"
        % (parameters[0], parameters[1], parameters[2], parameters[3], parameters[4]),
        fontsize=8,
    )

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")

    if show or not output_path:
        plt.show()
    else:
        plt.close(figure)


def main():
    args = parse_args()
    real_parts, imaginary_parts = load_coefficients(args.cms)
    parameters = load_parameters(args.parameters)
    plot_pattern(real_parts, imaginary_parts, parameters, output_path=args.output, show=args.show)


if __name__ == "__main__":
    main()
