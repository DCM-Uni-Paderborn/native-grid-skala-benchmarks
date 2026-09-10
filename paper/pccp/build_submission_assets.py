#!/usr/bin/env python3
"""Generate publication tables and figures from the frozen benchmark data."""
import csv
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
LC = REPO / 'benchmarks/Goldzak12'
X23 = REPO / 'benchmarks/X23-mini'
DIET = REPO / 'benchmarks/dietGMTKN55/production-p25/paper-common-70'
sys.path.insert(0, str(LC / 'scripts'))
from fit_eos import fit_curve, read_rows, write_rows
from fit_selected_eos import select_curves

LABELS = {'gapwxc-gth':'GX', 'gapw-ae':'AE', 'hybrid-direct':'HD', 'hybrid-one-center':'HOC'}
sources = []


def rows(path):
    sources.append(path)
    return read_rows(path)


def table(path, caption, spec, header, data, label=None, long=False):
    environment = 'longtable' if long else 'tabular'
    caption_text = '\\caption{'+caption+'}'
    if label:
        caption_text += '\\label{'+label+'}'
    text = ('{\\small\n\\begin{longtable}{'+spec+'}\n'+caption_text+'\\\\\n') if long else ('\\begin{table}[htbp]\n\\centering\\small\n'+caption_text+'\n')
    if not long:
        text += '\\begin{tabular}{'+spec+'}\n'
    text += '\\toprule\n' + ' & '.join(header) + '\\\\\n\\midrule\n'
    if long:
        text += '\\endfirsthead\n\\toprule\n' + ' & '.join(header) + '\\\\\n\\midrule\n\\endhead\n'
    text += '\n'.join(' & '.join(map(str,r)) + r' \\' for r in data)
    text += '\n\\bottomrule\n\\end{'+environment+'}\n' + ('}\n' if long else '\\end{table}\n')
    (HERE/path).write_text(text)


def molecular_tables():
    data = rows(DIET/'native-protocol-comparison.csv')
    assert len(data) == 70
    assert abs(np.mean([abs(float(r['gapwxc_gth_error_kcal_mol'])) for r in data])-3.5050583176051737)<1e-10
    table('molecular-reactions.tex', r'Complete common-set reaction energies in kcal mol$^{-1}$. Reference values are official benchmark energies. Protocol labels follow Sec.~\ref{sec:scope}. No reference-error outlier is omitted.', 'llrrrr',
          ['Subset','Reaction','Ref.','GX','HD','HOC'],
          [[r['subset'].replace('_',r'\_'),r['reaction_id']]+[f"{float(r[k]):.5f}" for k in ['reference_kcal_mol','gapwxc_gth_kcal_mol','hybrid_ae_gth_direct_kcal_mol','hybrid_ae_gth_one_center_kcal_mol']] for r in data],label='tab:molecular-reactions',long=True)
    cutoff = rows(REPO/'convergence/cutoff/aconf8-b200/results.csv')
    pairs = {}
    for row in cutoff:
        if row['converged'] != 'true':
            continue
        key = row['method'], int(float(row['cutoff_ry']))
        pairs.setdefault(key,{})[row['path'].split('-')[-1]] = float(row['total_energy_ha'])
    energies = {key:(p['gtg']-p['ttt'])*627.5094740631 for key,p in pairs.items() if set(p)=={'ttt','gtg'}}
    data = []
    for c in sorted({k[1] for k in energies}):
        if not all((m,c) in energies for m in ['GPW','GAPW_XC']):
            continue
        vals = [energies[m,c] for m in ['GPW','GAPW_XC']]
        errs = [energies[m,c]-energies[m,800] for m in ['GPW','GAPW_XC']]
        data.append([str(c)]+[f'{v:.6f}' for v in vals+errs])
    table('cutoff-table.tex','Conformational energy of the paired ACONF geometries in kcal mol$^{-1}$. Each error is relative to its own 800 Ry result. The relative cutoff is 60 Ry.', 'rrrrr',
          ['Cutoff/Ry','GPW','GAPW-XC',r'$\Delta$GPW',r'$\Delta$GAPW-XC'],data)
    grid = sorted(rows(REPO/'convergence/atom-grid/results/h2o-grid-results.csv'),
                  key=lambda r: (int(r['radial_grid']), int(r['lebedev_grid'])))
    reference = next(r for r in grid if r['radial_grid']=='200' and r['lebedev_grid']=='974')
    energy = float(reference['total_energy_ha'])
    data = [[r['radial_grid'],r['lebedev_grid'],f"{(float(r['total_energy_ha'])-energy)*1e6:+.3f}",f"{float(r['composite_electrons'])-8:+.3e}"] for r in grid if r['converged']=='true']
    table('atom-grid-table.tex','Self-consistent GPW water quadrature test at 800/60 Ry. Total-energy shifts are relative to 200/974 and are in microhartree. Electron-integral errors refer to eight explicit electrons.', 'rrrr', ['Radial','Angular',r'$\Delta E/\mu E_h$',r'$\Delta N/e$'],data)


def eos_observables():
    selection_path = LC/'protocol/paper-eos-selection.json'
    sources.append(selection_path)
    selection = json.loads(selection_path.read_text())
    curves = select_curves(rows(LC/'results/eos-selected-source.csv'),selection)
    base = {(r['method'],r['solid']):r for r in rows(LC/'results/eos-selected-fits.csv')}
    results = []
    cache = {}
    for key,(source,points) in sorted(curves.items()):
        b = base[key]
        cachekey = source,key[1]
        if cachekey not in cache:
            omitted = []
            for test,subset in [('omit-smallest',points[1:]),('omit-largest',points[:-1])]:
                pars,_,_ = fit_curve(subset)
                lo,hi = min(float(r['volume_A3']) for r in subset),max(float(r['volume_A3']) for r in subset)
                omitted.append({'test':test,'B0_prime':float(pars[3]),'delta_B0_prime':float(pars[3])-float(b['B0_prime']), 'minimum_bracketed':bool(lo<pars[1]<hi)})
            cache[cachekey]=omitted
        results.append({'method':key[0],'solid':key[1],'V0_A3_cell':float(b['V0_A3_cell']),
                        'B0_GPa':float(b['B0_GPa']),'B0_prime':float(b['B0_prime']),
                        'compressibility_GPa_inverse':1/float(b['B0_GPa']),
                        'max_endpoint_delta_B0_prime':max(abs(r['delta_B0_prime']) for r in cache[cachekey]),
                        'both_omitted_minima_bracketed':all(r['minimum_bracketed'] for r in cache[cachekey])})
    write_rows(HERE/'eos-additional-observables.csv',list(results[0]),results)
    (HERE/'eos-pressure-derivative-sensitivity.json').write_text(json.dumps({f'{a}/{b}':v for (a,b),v in cache.items()},indent=2)+'\n')
    table('eos-extra-table.tex',r'Additional predictions from the existing EOS data. $\kappa_0$ is in $10^{-3}$ GPa$^{-1}$. $V_0$ is the eight-atom cell volume in \AA$^3$. The endpoint column is the largest absolute change of $B_0^\prime$ after removing either endpoint. A dagger means at least one omitted-window minimum is no longer bracketed.',
          'llrrrrr',['Method','Solid',r'$V_0$',r'$B_0$/GPa',r'$B_0^\prime$',r'$10^3\kappa_0$',r'max. $|\Delta B_0^\prime|$'],
          [[LABELS[r['method']],r['solid'],f"{r['V0_A3_cell']:.4f}",f"{r['B0_GPa']:.2f}",f"{r['B0_prime']:.3f}",f"{1000*r['compressibility_GPa_inverse']:.3f}",f"{r['max_endpoint_delta_B0_prime']:.3f}"+(r'$^\dagger$' if not r['both_omitted_minima_bracketed'] else '')] for r in results],label='tab:eos-extra',long=True)
    print('Bprime range',min(r['B0_prime'] for r in results),max(r['B0_prime'] for r in results))
    print('Maximum Bprime endpoint sensitivity',max(results,key=lambda r:r['max_endpoint_delta_B0_prime']))


def copied_tables():
    for source in [LC/'paper/periodic-tables-si.tex',X23/'paper/x23-tables-si.tex']:
        sources.append(source)
        text = source.read_text().replace(r'\texttt{GAPW\_XC}', 'GAPW-XC')
        (HERE/source.name).write_text(text)


def eos_figure():
    data = rows(LC/'results/eos-literature-comparison-values.csv')
    selection_path = LC/'protocol/paper-eos-selection.json'
    sources.append(selection_path)
    selection = json.loads(selection_path.read_text())
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
    fig.savefig(HERE/'eos-literature-residuals.png', dpi=300)
    plt.close(fig)


def figures():
    plt.rcParams.update({'font.family':'sans-serif','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(figsize=(3.35,2.7))
    crystal_path=X23/'results/lattice-energies.json'
    sources.append(crystal_path)
    crystal_data=json.loads(crystal_path.read_text())
    assert crystal_data['accepted_base_cases']==24 and len(crystal_data['complete_pairs'])==12
    pairs={(r['method'],r['system']):r for r in crystal_data['complete_pairs']}
    methods=['gapwxc-gth','gapw-ae','gapw-gth-direct','gapw-gth-one-center']
    systems=['CO2','NH3','urea']
    reference=np.array([pairs[methods[0],s]['dmc_kjmol'] for s in systems])
    uncertainties=[pairs[methods[0],s]['dmc_statistical_uncertainty_kjmol'] for s in systems]
    values=[[pairs[m,s]['lattice_energy_kjmol'] for s in systems] for m in methods]
    for offset,(label,value,color,marker) in enumerate(zip(['GAPW-XC/GTH','GAPW-AE','GTH direct','GTH one centre'],values,['#007C83','#BD3F53','#566078','#D49B16'],['o','s','^','D'])):
        ax.scatter(np.arange(3)+(offset-1.5)*.12,np.asarray(value)-reference,s=30,label=label,color=color,marker=marker,zorder=3)
    ax.errorbar(np.arange(3),np.zeros(3),yerr=uncertainties,fmt='none',color='black',capsize=4,zorder=2)
    ax.axhline(0,color='#777777',lw=.7)
    ax.set_xticks(range(3),[r'CO$_2$',r'NH$_3$','Urea'])
    ax.set_ylabel(r'$E_{\rm latt}-E_{\rm DMC}$ / kJ mol$^{-1}$')
    ax.set_ylim(-18,3);ax.set_xlim(-.45,2.45)
    ax.legend(frameon=False,fontsize=7.5,loc='upper center',bbox_to_anchor=(.5,-.16),ncol=2)
    fig.tight_layout(pad=.4);fig.savefig(HERE/'crystal-lattice-errors.pdf');plt.close(fig)
    fig,ax=plt.subplots(figsize=(3.35,3.45));ax.set_xlim(0,10);ax.set_ylim(0,10);ax.axis('off')
    def box(x,y,text,width=2.8,color='#E7F2F2'):
        ax.text(x,y,text,ha='center',va='center',fontsize=8.5,bbox={'boxstyle':'square,pad=.5','fc':color,'ec':'#657575','lw':.7})
    box(2.5,8.8,'Smooth fields\nregular grids');box(7.5,8.8,'Hard minus soft\none-centre fields',color='#F7ECEB')
    box(5,6.4,r'Joint $\rho_\sigma,\,\nabla\rho_\sigma,\,\tau_\sigma$'+'\non atom-composite points')
    box(5,3.9,'Descriptors and nonlocal coupling\none reconstructed functional')
    box(5,1.4,r'$E_{\rm xc}$ and primitive-field derivatives'+'\nconsistent adjoint to CP2K')
    for start,end in [((2.5,7.9),(4,7)),((7.5,7.9),(6,7)),((5,5.6),(5,4.65)),((5,3.1),(5,2.15))]:
        ax.annotate('',xy=end,xytext=start,arrowprops={'arrowstyle':'->','color':'#007C83','lw':1.4})
    fig.tight_layout(pad=.05);fig.savefig(HERE/'native-reconstruction.pdf');plt.close(fig)


def toc_graphic():
    fig,ax=plt.subplots(figsize=(8/2.54,4/2.54))
    ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    ax.text(.5,.92,'Native Skala in CP2K',ha='center',va='center',fontsize=10,weight='bold')
    ax.text(.19,.60,r'$\widetilde{\rho}+\rho^h-\rho^s$'+'\n'+r'$\nabla\rho,\ \tau$',ha='center',va='center',fontsize=10)
    ax.text(.70,.60,'Joint functional\nGPW / GAPW',ha='center',va='center',fontsize=9,
            bbox={'boxstyle':'square,pad=.45','fc':'#E7F2F2','ec':'#007C83','lw':.8})
    ax.annotate('',xy=(.49,.60),xytext=(.36,.60),arrowprops={'arrowstyle':'->','color':'#007C83','lw':1.4})
    ax.plot([.07,.93],[.33,.33],color='#758287',lw=.6)
    ax.text(.27,.19,'Gas-phase\nmolecules',ha='center',va='center',fontsize=8.5)
    ax.text(.73,.19,'Condensed-phase\ncrystals',ha='center',va='center',fontsize=8.5)
    fig.subplots_adjust(left=.01,right=.99,bottom=.01,top=.99)
    fig.savefig(HERE/'toc-graphic.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
    plt.close(fig)


if __name__ == '__main__':
    molecular_tables();copied_tables();eos_observables();figures();eos_figure();toc_graphic()
    manifest={str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(sources))}
    (HERE/'source-data-sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')
