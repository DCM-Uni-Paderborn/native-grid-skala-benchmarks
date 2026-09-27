"""Figure routines used by the accepted PCCP manuscript."""
import math

PHASES = ("Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI",
          "XIII", "XIV", "XV", "XVII")
BASES = (("tzvpp", "TZVPP-ae", "#b34a32"),
         ("qzvpp", "QZVPP-ae", "#087e83"))

def ice_data(data):
    rows = data["rows"]
    assert tuple(row["phase"] for row in rows) == PHASES
    for row in rows:
        assert all(math.isfinite(row[key]) for key in
                   ("tzvpp", "qzvpp", "dmc_kjmol", "dmc_statistical_uncertainty_kjmol"))
        assert row["dmc_statistical_uncertainty_kjmol"] >= 0
    ih = rows[0]
    results = []
    for row in rows:
        result = {"phase": row["phase"]}
        for basis, _, _ in BASES:
            result[f"{basis}_absolute_error_kjmol"] = row[basis] - row["dmc_kjmol"]
            result[f"{basis}_relative_error_kjmol"] = (
                row[basis] - ih[basis] - (row["dmc_kjmol"] - ih["dmc_kjmol"]))
        results.append(result)
    stats = {kind: {basis: sum(abs(row[f"{basis}_{kind}_error_kjmol"])
                              for row in selected) / len(selected)
                    for basis, _, _ in BASES}
             for kind, selected in (("absolute", results), ("relative", results[1:]))}
    return results, stats

def ice_figure(results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator

    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 9,
                         "pdf.fonttype": 42, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.linewidth": 0.6}):
        fig, ax = plt.subplots(figsize=(6.8, 2.9))
        fig.subplots_adjust(left=0.14, right=0.985, bottom=0.19, top=0.84)
        selected = results[1:]
        positions = list(range(len(selected)))
        for offset, (basis, label, color) in zip((-0.18, 0.18), BASES, strict=True):
            values = [row[f"{basis}_relative_error_kjmol"] for row in selected]
            ax.bar([x + offset for x in positions], values, width=0.32,
                   color=color, label=label, zorder=3)
        ax.axhline(0, color="#444444", linewidth=0.8, zorder=2)
        ax.set_xticks(positions, [row["phase"] for row in selected])
        ax.set_xlim(-0.65, len(selected) - 0.35)
        ax.set_ylim(-9.7, 1.4)
        ax.yaxis.set_major_locator(MultipleLocator(2))
        ax.grid(axis="y", color="#dddddd", linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)
        ax.set_ylabel("Ih-relative energy error\n(kJ mol$^{-1}$)", labelpad=6)
        ax.set_xlabel("Ice phase", labelpad=7)
        ax.tick_params(axis="both", length=3, width=0.6)
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.985, 0.997),
                   ncol=2, frameon=False, handlelength=1.3, columnspacing=1.5)
        return fig, ax


def molecular_figure(results, counts):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator

    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 8,
                         "pdf.fonttype": 42, "axes.spines.top": False,
                         "axes.spines.right": False}):
        fig, axes = plt.subplots(1, 2, figsize=(6.5, 8.5), sharey=True)
        fig.subplots_adjust(left=0.185, right=0.985, bottom=0.07, top=0.93, wspace=0.16)
        all_values = [row[f"improvement_{name}_kj_mol"]
                      for row in results for name in ("reference", "gauxc")]
        low = 5 * math.floor(min(all_values) / 5)
        high = 5 * math.ceil(max(all_values) / 5)
        positions = list(range(len(results)))
        names = [row["reaction"] for row in results]
        for ax, name, title in zip(axes, ("reference", "gauxc"),
                                   ("(a) Benchmark reference", "(b) GauXC-AE"), strict=True):
            values = [row[f"improvement_{name}_kj_mol"] for row in results]
            colors = ["#087e83" if value > 0 else "#b34a32" for value in values]
            ax.barh(positions, values, color=colors, height=0.68, zorder=3)
            ax.scatter(values, positions, c=colors, s=5, zorder=4)
            ax.axvline(0, color="#555555", lw=0.7, zorder=2)
            ax.set_xlim(low, high)
            ax.set_ylim(len(results) - 0.35, -0.8)
            ax.xaxis.set_major_locator(MultipleLocator(10))
            ax.xaxis.set_minor_locator(MultipleLocator(5))
            ax.grid(axis="x", color="#dedede", lw=0.45)
            ax.set_axisbelow(True)
            ax.set_title(f"{title}\n{counts[name]}/62 smaller deviations", fontsize=9, pad=9)
            ax.set_xlabel("Decrease in absolute deviation\n(kJ mol$^{-1}$)", fontsize=8)
            ax.tick_params(axis="y", length=0, pad=4, labelsize=7)
            ax.tick_params(axis="x", labelsize=8)
        axes[0].set_yticks(positions, names)
        axes[1].tick_params(axis="y", labelleft=False)
        return fig, axes


def eos_figure(data, selection):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    solids, methods = selection['solids'], selection['methods']
    refs = {r['solid']: r for r in data if r['source'] == 'Goldzak2022' and r['method'] == 'experiment'}
    fig, axes = plt.subplots(1, 2, figsize=(7, 3.8), sharey=True)
    palette = ['#007C83', '#BD3F53', '#526936', '#80619E']
    labels = {'gapwxc-gth': 'GAPW-XC/GTH', 'gapw-ae': 'GAPW-AE',
              'hybrid-direct': 'Mixed direct', 'hybrid-one-center': 'Mixed one centre'}
    for ax, prop, label in zip(axes, ['a_A', 'B0_GPa'],
                                [r'$a_0-a_{\rm ref}$ / $\AA$', r'$B_0-B_{\rm ref}$ / GPa']):
        for i, method in enumerate(methods):
            values = {r['solid']: float(r[prop]) for r in data if r['source'] == 'NativeSkala' and r['method'] == method}
            ax.scatter([values[s]-float(refs[s][prop]) for s in solids],
                       np.arange(len(solids))+(i-1.5)*.12, s=18,
                       label=labels[method], color=palette[i], marker=['o','s','^','x'][i], zorder=3)
        scs = {r['solid']: float(r[prop]) for r in data if r['source'] == 'Goldzak2022' and r['method'] == 'SCS-MP2'}
        ax.scatter([scs[s]-float(refs[s][prop]) for s in solids], np.arange(len(solids)),
                   label='SCS-MP2', marker='D', s=13, color='#303030', zorder=3)
        ax.axvline(0, color='#777777', lw=.7)
        ax.grid(axis='x', color='#E2E2E2', lw=.5)
        ax.set_xlabel(label)
        ax.tick_params(labelsize=8)
        ax.set_yticks(np.arange(len(solids)), solids)
    axes[0].invert_yaxis()
    axes[0].set_xticks([-.15, -.10, -.05, 0, .025])
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=3, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0,.13,1,1), pad=.6)
    return fig


def reconstruction_figure():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(3.35, 3.45))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')
    def box(x, y, text, color='#E7F2F2'):
        ax.text(x, y, text, ha='center', va='center', fontsize=8.5,
                bbox={'boxstyle': 'square,pad=.5', 'fc': color, 'ec': '#657575', 'lw': .7})
    box(2.5, 8.8, 'Smooth fields\nregular grids')
    box(7.5, 8.8, 'Hard minus soft\none-centre fields', '#F7ECEB')
    box(5, 6.4, r'Joint $\rho_\sigma,\,\nabla\rho_\sigma,\,\tau_\sigma$'+'\non atom-composite points')
    box(5, 3.9, 'Descriptors and nonlocal coupling\none reconstructed functional')
    box(5, 1.4, r'$E_{\rm xc}$ and primitive-field derivatives'+'\nconsistent adjoint to CP2K')
    for start, end in [((2.5,7.9),(4,7)), ((7.5,7.9),(6,7)), ((5,5.6),(5,4.65)), ((5,3.1),(5,2.15))]:
        ax.annotate('', xy=end, xytext=start, arrowprops={'arrowstyle':'->','color':'#007C83','lw':1.4})
    fig.tight_layout(pad=.05)
    return fig
