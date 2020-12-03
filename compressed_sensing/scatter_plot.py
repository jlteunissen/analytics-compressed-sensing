
from bokeh.models import TapTool, CustomJS, ColumnDataSource, HoverTool, ColumnDataSource
from bokeh.io import show, output_notebook
from bokeh.plotting import figure, show
from bokeh.embed import components
import itertools
from bokeh.palettes import Dark2_5 as palette

output_notebook()

def show_scatter_plot(xs, ys, data_point_labels=None, x_label=None, y_label=None, legend=None, unit=None):
    # if xs ist not list of lists/arrays make it so, as later the function iterates over xs and ys
    if not isinstance(xs[0], (list, np.ndarray)):
        xs = [xs]
        ys = [ys]
    # make sure that xs and ys ist list (of lists/arrays) as later the function will
    # do the list operation xs+ys
    elif not isinstance(xs, list) or not isinstance(ys, list):
        xs = list(xs)
        ys = list(ys)

    if unit is None:
        unit = ''

    hover = HoverTool(
        tooltips="""
            <div>
                <div>
                    <span style="font-size: 15px; font-weight: bold;">@data_point_labels</span>
                </div>
                <div >
                    <span style="font-size: 10px;">Abs. error = @abs_error %s</span><br>
                </div>
                <div>
                    <span style="font-size: 10px;">Location:</span>
                    <span style="font-size:  10px; color: #696;">($x, $y)</span>
                </div>
            </div>
            """ % unit
    )

    colors = itertools.cycle(palette)

    p = figure(plot_width=600, plot_height=300, tools=[hover, "box_zoom", "pan", "reset"],
               x_axis_label=x_label, y_axis_label=y_label)

    # plot reference diagonal
    xy_min = min([min(arr) for arr in xs + ys])
    xy_max = max([max(arr) for arr in xs + ys])
    p.line([xy_min, xy_max], [xy_min, xy_max])

    for i, color in zip(range(len(xs)), colors):
        source = ColumnDataSource(
            data=dict(
                x=xs[i],
                y=ys[i],
                data_point_labels=data_point_labels[i],
                abs_error=abs(np.array(xs[i]) - np.array(ys[i]))
            )
        )

        p.circle('x', 'y', size=8, source=source, legend=legend[i], color=color)
    p.legend.location = 'top_left'
    show(p)
