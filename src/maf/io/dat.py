"""Plain numeric DAT tables with quantity/unit headers shared by plots."""
from pathlib import Path
import json
import re
import numpy as np


def observable_axes(definition, units):
    kind=definition['kind']
    coordinate={'dndy':'y','net_proton':'y','dndeta':'eta','pt':'p_T','invariant_pt':'p_T',
                'mt':'m_T','invariant_mt':'m_T','mt_minus_m0':'m_T-m_0'}.get(kind,'index')
    if kind in ('v1','v2','vn'):
        coordinate={'y':'y','eta':'eta','pt':'p_T'}[definition.get('axis','y')]
    quantity={'yield':'<N>','charged_multiplicity':'<N_ch>','dndy':'dN/dy','dndeta':'dN_ch/deta',
              'net_proton':'d(N_p-N_pbar)/dy','pt':'dN/dp_T','mt':'dN/dm_T',
              'mt_minus_m0':'dN/d(m_T-m_0)','invariant_pt':'(1/(2*pi*p_T)) d^2N/(dp_T dy)',
              'invariant_mt':'(1/(2*pi*m_T)) d^2N/(dm_T dy)','mean_pt':'<p_T>',
              'mean_mt_minus_m0':'<m_T-m_0>','ratio':'N_numerator/N_denominator',
              'v1':'v_1','v2':'v_2'}.get(kind,kind)
    if kind in ('v1','v2','vn'):
        harmonic=definition.get('harmonic') or (int(kind[1:]) if kind != 'vn' else None)
        method=(definition.get('flow') or {}).get('method')
        suffix={'reaction_plane':'RP','event_plane':'EP','cumulant_2':'2','cumulant_4':'4'}.get(method)
        quantity=f'v_{harmonic}' + (f'{{{suffix}}}' if suffix else '')
    density=definition.get('density')
    if density:
        axis={'y':'y','eta':'eta'}[density]
        if kind in ('yield','charged_multiplicity'):
            quantity=('dN_ch/d' if kind=='charged_multiplicity' else 'dN/d')+axis
        elif kind in ('pt','mt','mt_minus_m0'):
            quantity=f'd^2N/(d({coordinate}) d{axis})'
    return {'x':{'quantity':coordinate,'unit':'GeV' if coordinate in ('p_T','m_T','m_T-m_0') else '1'},
            'y':{'quantity':quantity,'unit':units}}


def _math_quantity(quantity):
    """Translate canonical DAT quantities to MathText; preserve custom prose."""
    coordinates={'y':'y', 'eta':r'\eta', 'p_T':r'p_{T}',
                 'm_T':r'm_{T}', 'm_T-m_0':r'm_{T}-m_{0}'}
    quantities=dict(coordinates)
    quantities.update({
        '<N>':r'\langle N\rangle', '<N_ch>':r'\langle N_{\mathrm{ch}}\rangle',
        '<p_T>':r'\langle p_{T}\rangle', '<m_T-m_0>':r'\langle m_{T}-m_{0}\rangle',
        'd(N_p-N_pbar)/dy':r'\frac{\mathrm{d}(N_{p}-N_{\bar{p}})}{\mathrm{d}y}',
        'N_numerator/N_denominator':r'\frac{N_{\mathrm{numerator}}}{N_{\mathrm{denominator}}}',
    })
    for axis in ('y','eta'):
        for count,symbol in (('N','N'),('N_ch',r'N_{\mathrm{ch}}')):
            quantities[f'd{count}/d{axis}']=rf'\frac{{\mathrm{{d}}{symbol}}}{{\mathrm{{d}}{coordinates[axis]}}}'
    for coordinate in ('p_T','m_T','m_T-m_0'):
        symbol=coordinates[coordinate]
        differential=rf'\mathrm{{d}}({symbol})' if coordinate=='m_T-m_0' else rf'\mathrm{{d}}{symbol}'
        quantities[f'dN/d{coordinate}']=rf'\frac{{\mathrm{{d}}N}}{{{differential}}}'
        quantities[f'dN/d({coordinate})']=rf'\frac{{\mathrm{{d}}N}}{{{differential}}}'
        for axis in ('y','eta'):
            quantities[f'd^2N/(d({coordinate}) d{axis})']=rf'\frac{{\mathrm{{d}}^{{2}}N}}{{\mathrm{{d}}({symbol})\,\mathrm{{d}}{coordinates[axis]}}}'
        if coordinate!='m_T-m_0':
            key=f'(1/(2*pi*{coordinate})) d^2N/(d{coordinate} dy)'
            quantities[key]=rf'\frac{{1}}{{2\pi {symbol}}}\frac{{\mathrm{{d}}^{{2}}N}}{{{differential}\,\mathrm{{d}}y}}'
    if quantity in quantities:
        return quantities[quantity]
    flow=re.fullmatch(r'v_(\d+)(?:\{(RP|EP|2|4)\})?',quantity)
    if flow:
        label=rf'v_{{{flow[1]}}}'
        if flow[2]:
            method=rf'\mathrm{{{flow[2]}}}' if flow[2].isalpha() else flow[2]
            label+=rf'\{{{method}\}}'
        return label
    subscript=re.fullmatch(r'([A-Za-z])_(?:\{([A-Za-z0-9]+)\}|([A-Za-z0-9]+))',quantity)
    if subscript:
        index=subscript[2] or subscript[3]
        index=rf'\mathrm{{{index}}}' if index.isalpha() and len(index)>1 else index
        return rf'{subscript[1]}_{{{index}}}'
    return None


def axis_label(axis):
    """Render plot labels with MathText, retaining plain quantities in DAT."""
    math=_math_quantity(axis['quantity'])
    quantity=f'${math}$' if math is not None else axis['quantity']
    if axis['unit']=='1':
        return quantity
    unit=axis['unit']
    power=re.fullmatch(r'([A-Za-z]+)\^([+-]?\d+)',unit)
    if power:
        unit=rf'$\mathrm{{{power[1]}}}^{{{power[2]}}}$'
    return f'{quantity} [{unit}]'


def write_dat(path, name, data, *, signature='', run_pattern=None):
    axes=data.get('axes') or observable_axes(data['definition'],data['units'])
    values=np.asarray(data['value'],float)
    errors=np.asarray(data['stat_error'],float)
    bins=data['bins']
    if bins is None:
        centers=np.arange(len(values),dtype=float)
        low=high=np.full(len(values),np.nan)
    else:
        edges=np.asarray(bins,float)
        low,high=edges[:-1],edges[1:]
        centers=(low+high)/2
    x,y=axes['x'],axes['y']
    header=[f'observable: {name}',f'x_quantity: {x["quantity"]}',f'x_unit: {x["unit"]}',
            f'y_quantity: {y["quantity"]}',f'y_unit: {y["unit"]}',
            f'analysis_signature: {signature}',
            'definition: '+json.dumps(data['definition'],sort_keys=True,allow_nan=False),
            'columns: '+ '\t'.join([f'{x["quantity"]}[{x["unit"]}]',f'{x["quantity"]}_low[{x["unit"]}]',
                f'{x["quantity"]}_high[{x["unit"]}]',f'{y["quantity"]}[{y["unit"]}]',f'stat_error[{y["unit"]}]'])]
    if run_pattern is not None:
        header.append('run_pattern: '+json.dumps(run_pattern,ensure_ascii=False))
    np.savetxt(path,np.column_stack((centers,low,high,values,errors)),fmt='%.17g',delimiter='\t',header='\n'.join(header))


def read_dat(path):
    path=Path(path)
    header={}
    with path.open() as stream:
        for line in stream:
            if not line.startswith('#'):
                break
            key,separator,value=line[1:].strip().partition(':')
            if separator:
                header[key]=value.strip()
    required={'observable','x_quantity','x_unit','y_quantity','y_unit','definition'}
    if not required.issubset(header) or any(not header[k] for k in required):
        raise ValueError(f'Missing DAT quantity/unit/definition header: {path}')
    definition=json.loads(header['definition'])
    table=np.loadtxt(path,ndmin=2)
    if table.shape[1]!=5 or len(table)==0:
        raise ValueError(f'DAT requires center, lower edge, upper edge, value and statistical error: {path}')
    if definition.get('bins') is not None:
        bins=np.asarray(definition['bins'],float)
        if len(bins)!=len(table)+1 or not np.allclose(table[:,1],bins[:-1]) or not np.allclose(table[:,2],bins[1:]):
            raise ValueError(f'DAT edges do not match the saved definition: {path}')
    if np.any(table[:,4]<0):
        raise ValueError(f'Negative DAT statistical error: {path}')
    data={'definition':definition,'bins':definition.get('bins'),'x':table[:,0],
          'value':table[:,3],'stat_error':table[:,4],'units':header['y_unit'],
          'axes':{'x':{'quantity':header['x_quantity'],'unit':header['x_unit']},
                  'y':{'quantity':header['y_quantity'],'unit':header['y_unit']}},
          'analysis_signature':header.get('analysis_signature','')}
    if 'run_pattern' in header:
        data['run_pattern']=json.loads(header['run_pattern'])
    return header['observable'],data


def load_summary_dat(path):
    """DAT values and header labels take precedence for saved-result plotting."""
    path=Path(path)
    summary=json.loads(path.read_text())
    for name,data in summary['observables'].items():
        dat=path.parent/f'{name}.dat'
        if not dat.exists():
            continue  # Legacy JSON/CSV results remain plottable.
        saved_name,saved=read_dat(dat)
        if saved_name!=name or (saved['analysis_signature'] and saved['analysis_signature']!=summary['analysis_signature']):
            raise ValueError(f'DAT does not belong to this result: {dat}')
        data.update(saved)
    return summary
