import numpy as np
import matplotlib.pyplot as plt
import os


def compute_fix_ratio(flt_file_path):
    """
    Extract the AmbStatus field from a .flt file and calculate the fixing rate (Fixed / Total).
    """
    total = 0
    fixed = 0

    with open(flt_file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("#") or len(line.strip()) == 0:
                continue
            parts = line.strip().split()
            if len(parts) < 17:
                continue
            status = parts[16]
            if status in ["Fixed", "Float"]:
                total += 1
                if status == "Fixed":
                    fixed += 1

    if total == 0:
        return None
    return fixed / total * 100  # Return fixing rate as percentage


def plot_pseudorange_data(filename, prefix="G", save_path=None, fix_ratio=None):
    """
    Read STEC data and generate a plot. Optionally display fixing rate if available.
    """
    data = []

    with open(filename, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue

        parts = line.split()
        values = parts[2:] if parts[0] == "Epoch" else parts

        row_data = []
        for x in values:
            try:
                val = float(x)
                row_data.append(np.nan if val == 0.0 else val)
            except ValueError:
                row_data.append(np.nan)
        data.append(row_data)

    try:
        data = np.array(data, dtype=float)
    except:
        return

    if data.shape[0] == 0 or data.shape[1] == 0:
        return

    x = np.arange(data.shape[0]) / 120.0  # Time in hours
    num_satellites = data.shape[1]
    colors = plt.get_cmap('rainbow')(np.linspace(0, 1, num_satellites))

    plt.figure(figsize=(12, 6))
    label_prefix = "E" if prefix == "EX" else prefix

    for i in range(num_satellites):
        prn = f"{label_prefix}{i + 1:02d}"
        plt.plot(x, data[:, i], color=colors[i], label=prn)

    plt.xlabel("Time (hours)")
    plt.ylabel("STEC (TECU)")
    plt.xlim(0, 24)
    plt.xticks(np.arange(0, 25, 2))
    plt.title(os.path.splitext(os.path.basename(filename))[0])
    plt.grid(True)
    plt.legend(loc='upper left', bbox_to_anchor=(1.01, 1), ncol=2, fontsize=8)
    plt.tight_layout()

    if fix_ratio is not None:
        x_min, x_max = plt.xlim()
        y_min, y_max = plt.ylim()
        plt.annotate(
            f"Fixing Rate: {fix_ratio:.2f}%",
            xy=(x_max, y_min),
            xytext=(x_max - 0.1, y_min - (y_max - y_min) * 0.06),
            textcoords='data',
            fontsize=10,
            ha='right',
            va='top'
        )

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300)
        print(f"[✓] Image saved: {save_path}")

    plt.close()


def batch_plot(site_list, input_path, output_path, year, doy):
    """
    Batch process and plot STEC data for multiple stations and GNSS systems.
    """
    systems = {
        "GPS": "G",
        "GLO": "R",
        "GAL": "E",
        "GALX": "EX",
        "BDS": "C"
    }

    for site in site_list:
        # Try to locate .flt file containing this station name
        fix_ratio = None
        for file in os.listdir(input_path):
            if file.endswith(".flt") and site in file:
                flt_path = os.path.join(input_path, file)
                fix_ratio = compute_fix_ratio(flt_path)
                break

        for suffix, prefix in systems.items():
            txt_path = os.path.join(input_path, f"{site}_{suffix}.txt")
            subfolder = os.path.join(output_path, "ppp_stec_plot", f"{year}_{doy:03d}", site)
            png_path = os.path.join(subfolder, f"{site}_{suffix}.png")

            if os.path.exists(txt_path):
                plot_pseudorange_data(txt_path, prefix=prefix, save_path=png_path, fix_ratio=fix_ratio)


if __name__ == "__main__":
    site_list = ["HKWS", "ABPO"]
    input_path = "../data_ppp"
    output_path = "../output"
    year = 2023
    doy = 330

    batch_plot(site_list, input_path, output_path, year, doy)
