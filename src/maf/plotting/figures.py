"""Optional Matplotlib layer. Numerical analysis never imports Matplotlib."""
from pathlib import Path
import numpy as np
from ..io.dat import observable_axes, axis_label


def _plt():
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError("Install maf[plot]") from exc
    return plt


def _simulation_label(summary, observable):
    return (summary.get("run_pattern") or summary.get("simulation_label") or
            (str(summary["group"]) if summary.get("group") else observable))


def distribution(results, observable, *, labels=None, logy=False):
    if not results:
        raise ValueError("No results to plot")
    summaries=[result.summary() if hasattr(result,"summary") else result for result in results]
    first=summaries[0]["observables"][observable]
    axes=first.get("axes") or observable_axes(first["definition"],first["units"])
    for summary in summaries[1:]:
        current=summary["observables"][observable]
        current_axes=current.get("axes") or observable_axes(current["definition"],current["units"])
        if (current["units"],current["definition"]["kind"],current["definition"].get("axis","y")) != (first["units"],first["definition"]["kind"],first["definition"].get("axis","y")):
            raise ValueError("Cannot overlay different normalizations or axes; select compatible current results")
        if current_axes!=axes:
            raise ValueError("Cannot overlay different quantities or axis units")
    plt=_plt()
    fig,ax=plt.subplots(layout="constrained")
    for i,summary in enumerate(summaries):
        data=summary["observables"][observable]
        bins=data["bins"]
        x=data.get("x")
        if x is None:
            x=(np.asarray(bins[:-1])+np.asarray(bins[1:]))/2 if bins else np.arange(len(data["value"]))
        y=np.asarray(data["value"],float)
        error=np.asarray(data["stat_error"],float)
        ax.errorbar(x,y,yerr=error,marker="o",linewidth=1,label=labels[i] if labels else _simulation_label(summary,observable))
    ax.set(xlabel=axis_label(axes["x"]),ylabel=axis_label(axes["y"]))
    if logy:
        ax.set_yscale("log")
    ax.legend()
    return fig,ax


def all_distributions(summary,output,*,format="png",logy=False):
    """One saved graph per computed observable, closing figures as we go."""
    if format not in ("png","pdf"):
        raise ValueError("Batch plot format must be png or pdf")
    output=Path(output)
    output.mkdir(parents=True,exist_ok=True)
    plt=_plt()
    spectra={"pt","mt","mt_minus_m0","invariant_pt","invariant_mt"}
    written=[]
    for name,data in summary["observables"].items():
        if Path(name).name!=name or name in (".",".."):
            raise ValueError("Observable name must be a safe filename component")
        automatic_log=data["definition"]["kind"] in spectra and np.any(np.asarray(data["value"],float)>0)
        fig,_=distribution([summary],name,logy=logy or automatic_log)
        try:
            path=output/f"{name}.{format}"
            fig.savefig(path,metadata={"Creator":"maf"})
            written.append(path)
        finally:
            plt.close(fig)
    return written


def comparison(result,observable,dataset):
    from ..experiment.comparison import compare
    values=compare(result,observable,dataset)
    plt=_plt()
    fig,(ax,ratio)=plt.subplots(2,1,sharex=True,height_ratios=[3,1],layout="constrained")
    x=(dataset.low+dataset.high)/2
    xe=np.vstack((x-dataset.low,dataset.high-x))
    ax.errorbar(x,values["model"],yerr=values["model_stat"],label=result.run_pattern or "Model",marker="o")
    ax.errorbar(x,dataset.value,yerr=dataset.stat.T if dataset.stat is not None else None,xerr=xe,fmt="s",label=dataset.metadata["experiment"])
    if dataset.systematic is not None:
        ax.fill_between(x,dataset.value-dataset.systematic[:,0],dataset.value+dataset.systematic[:,1],alpha=.2,label="Data systematic")
    ratio.errorbar(x,values["ratio"],yerr=values["ratio_model_stat"],fmt="o",label="Model stat.")
    valid=dataset.value!=0
    if dataset.stat is not None:
        ratio.errorbar(x[valid],np.ones(np.count_nonzero(valid)),yerr=(dataset.stat[valid]/np.abs(dataset.value[valid,None])).T,fmt="none",color="black")
    ratio.axhline(1,color="grey",linewidth=.8)
    axes=observable_axes(dataset.metadata["definition"],dataset.metadata["definition"]["units"])
    ratio.set(xlabel=axis_label(axes["x"]),ylabel="Model/data")
    ax.set_ylabel(axis_label(axes["y"]))
    ax.legend()
    return fig,(ax,ratio)
