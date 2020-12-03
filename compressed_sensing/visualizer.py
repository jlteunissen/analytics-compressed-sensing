import plotly.graph_objects as go
import ipywidgets as widgets
from jupyter_jsmol import JsmolView
import numpy as np
from IPython.display import display, HTML, FileLink


class Visualizer:

    def __init__(self, df_D, sisso, D_selected_df):

        self.df_D = df_D
        self.sisso = sisso
        self.D_selected_df = D_selected_df
        self.marker_size = 7
        self.marker_symbol = 'circle'
        self.current_features = [0, 1]
        self.font_size = 12
        self.cross_size = 15
        self.line_width = 1
        self.font_families = ['Source Sans Pro',
                              'Helvetica',
                              'Open Sans',
                              'Times New Roman',
                              'Arial',
                              'Verdana',
                              'Courier New',
                              'Comic Sans MS',
                              ]
        self.line_styles = ["dash",
                            "solid",
                            "dot",
                            "longdash",
                            "dashdot",
                            "longdashdot"]
        self.bg_color = 'rgba(229,236,246, 0.5)'
        self.zb_color = "#EB8273"
        self.rs_color = "rgb(138, 147, 248)"
        self.total_features = sisso.n_nonzero_coefs
        self.features = []
        for i in range(self.total_features):
            self.features.append(df_D.columns[sisso.l0_selected_indices[self.total_features - 1]][i])
        self.coefficients = []
        for i in range(self.total_features):
            self.coefficients.append(sisso.coefs[sisso.l0_selected_indices[self.total_features - 1][i]])
        self.intercept = sisso.intercept
        self.line_x = np.linspace(D_selected_df[self.features[0]].min(), D_selected_df[self.features[0]].max(), 1000)
        self.line_y = self.f_x(self.line_x)
        self.fig = go.FigureWidget()
        self.viewer_l = JsmolView()
        self.viewer_r = JsmolView()
        self.bg_toggle = True
        self.compounds_list = df_D.index.tolist()

        self.text_RS = []
        for material in D_selected_df['Chem Formula'].tolist():
            self.text_RS.append(material + ' - RS structure')
        self.text_ZB = []
        for material in D_selected_df['Chem Formula'].tolist():
            self.text_ZB.append(material + ' - ZB structure')

        custom_RS = np.dstack((D_selected_df.loc[D_selected_df['Structure'] == 'RS']['energy_diff'],
                               D_selected_df.loc[D_selected_df['Structure'] == 'RS']['P_predict']))[0]
        custom_ZB = np.dstack((D_selected_df.loc[D_selected_df['Structure'] == 'ZB']['energy_diff'],
                               D_selected_df.loc[D_selected_df['Structure'] == 'ZB']['P_predict']))[0]

        # the final plot is the sum of two traces, respectively containing the RS vs ZB materials
        x_RS = D_selected_df.loc[D_selected_df['Structure'] == 'RS'][self.features[0]].to_numpy()
        y_RS = D_selected_df.loc[D_selected_df['Structure'] == 'RS'][self.features[1]].to_numpy()
        x_ZB = D_selected_df.loc[D_selected_df['Structure'] == 'ZB'][self.features[0]].to_numpy()
        y_ZB = D_selected_df.loc[D_selected_df['Structure'] == 'ZB'][self.features[1]].to_numpy()
        self.fig.add_trace(
            (
                go.Scatter(
                    mode='markers',
                    x=x_RS,
                    y=y_RS,
                    customdata=custom_RS,
                    text=D_selected_df.loc[D_selected_df['Structure'] == 'RS'][['Chem Formula']],
                    hovertemplate=
                    r"<b>%{text}</b><br><br>" +
                    "x axis: %{x:,.2f}<br>" +
                    "y axis: %{y:,.2f}<br>" +
                    "ΔE reference:  %{customdata[0]:,.4f}<br>" +
                    "ΔE predicted:  %{customdata[1]:,.4f}<br>",
                    name='RS',
                    marker=dict(color=self.rs_color)
                )
            ))
        self.fig.add_trace(
            (
                go.Scatter(
                    mode='markers',
                    x=x_ZB,
                    y=y_ZB,
                    customdata=custom_ZB,
                    text=D_selected_df.loc[D_selected_df['Structure'] == 'ZB'][['Chem Formula']],
                    hovertemplate=
                    r"<b>%{text}</b><br><br>" +
                    "x axis: %{x:,.2f}<br>" +
                    "y axis: %{y:,.2f}<br>" +
                    "ΔE reference:  %{customdata[0]:,.4f}<br>" +
                    "ΔE predicted:  %{customdata[1]:,.4f}<br>",
                    name='ZB',
                    marker=dict(color=self.zb_color)
                )
            ))
        self.fig.add_trace(
            (
                go.Scatter(
                    x=self.line_x,
                    y=self.line_y,
                    line=dict(color='Grey', width=1, dash=self.line_styles[0]),
                    name=r'Classification' + '<br>' + 'line',
                )
            )
        )

        x_min = min(min(x_RS), min(x_ZB))
        y_min = min(min(y_RS), min(y_ZB))
        x_max = max(max(x_RS), max(x_ZB))
        y_max = max(max(y_RS), max(y_ZB))
        x_delta = 0.05 * abs(x_max - x_min)
        y_delta = 0.05 * abs(y_max - y_min)
        self.fig.update_layout(
            plot_bgcolor=self.bg_color,
            font=dict(
                size=int(self.font_size),
                family=self.font_families[0]
            ),
            xaxis_title=self.features[0],
            yaxis_title=self.features[1],
            xaxis_range=[x_min - x_delta, x_max + x_delta],
            yaxis_range=[y_min - y_delta, y_max + y_delta],
            hoverlabel=dict(
                bgcolor="white",
                font_size=16,
                font_family="Rockwell"
            ),
            width=800,
            height=400,
            margin=dict(
                l=50,
                r=50,
                b=70,
                t=20,
                pad=4
            ),
        )

        self.scatter_RS = self.fig.data[0]
        self.scatter_ZB = self.fig.data[1]
        self.scatter_line = self.fig.data[2]
        self.RS_npoints = len(D_selected_df.loc[D_selected_df['Structure'] == 'RS'])
        self.ZB_npoints = len(D_selected_df.loc[D_selected_df['Structure'] == 'ZB'])
        self.RS_symbols = [self.marker_symbol] * self.RS_npoints
        self.ZB_symbols = [self.marker_symbol] * self.ZB_npoints
        self.RS_sizes = [self.marker_size] * self.RS_npoints
        self.ZB_sizes = [self.marker_size] * self.ZB_npoints
        self.update_markers()

        self.widg_featx = widgets.Dropdown(
            description='x-axis',
            options=self.features,
            value=self.features[0]
        )
        self.widg_featy = widgets.Dropdown(
            description='y-axis',
            options=self.features,
            value=self.features[1]
        )
        self.widg_featmarker = widgets.Dropdown(
            description="Marker",
            options=['Default size'] + self.features,
            value='Default size',
        )
        self.widg_compound_text_l = widgets.Combobox(
            placeholder='...',
            description='Compound:',
            options=self.compounds_list,
            disabled=False,
            layout=widgets.Layout(width='200px')
        )
        self.widg_compound_text_r = widgets.Combobox(
            placeholder='...',
            description='Compound:',
            options=self.compounds_list,
            disabled=False,
            layout=widgets.Layout(width='200px')
        )
        self.widg_display_button_l = widgets.Button(
            description="Display",
            layout=widgets.Layout(width='100px')
        )
        self.widg_display_button_r = widgets.Button(
            description="Display",
            layout=widgets.Layout(width='100px')
        )
        self.widg_checkbox_l = widgets.Checkbox(
            value=True,
            indent=False,
            layout=widgets.Layout(width='50px')
        )
        self.widg_checkbox_r = widgets.Checkbox(
            value=False,
            indent=False,
            layout=widgets.Layout(width='50px'),
        )
        self.widg_markersize = widgets.BoundedIntText(
            placeholder=str(self.marker_size),
            description='Marker size',
            value=str(self.marker_size)
        )
        self.widg_crosssize = widgets.BoundedIntText(
            placeholder=str(self.cross_size),
            description='Cross size',
            value=str(self.cross_size)
        )
        self.widg_fontsize = widgets.BoundedIntText(
            placeholder=str(self.font_size),
            description='Font size',
            value=str(self.font_size)
        )
        self.widg_linewidth = widgets.BoundedIntText(
            placeholder=str(self.line_width),
            description='Line width',
            value=str(self.line_width)
        )
        self.widg_linestyle = widgets.Dropdown(
            options=self.line_styles,
            description='Line style',
            value=self.line_styles[0],
        )
        self.widg_fontfamily = widgets.Dropdown(
            options=self.font_families,
            description='Font family',
            value=self.font_families[0]
        )
        self.widg_bgtoggle_button = widgets.Button(
            description='Toggle on/off background',
            layout=widgets.Layout(width='300px'),
        )
        self.widg_bgcolor = widgets.Text(
            placeholder=str(self.bg_color),
            description='Color',
            value=str(self.bg_color),
        )
        self.widg_rscolor = widgets.Text(
            placeholder=str(self.rs_color),
            description='RS color',
            value=str(self.rs_color),
        )
        self.widg_zbcolor = widgets.Text(
            placeholder=str(self.zb_color),
            description='ZB color',
            value=str(self.zb_color),
        )
        self.widg_markersymbol = widgets.Text(
            placeholder=str(self.marker_symbol),
            description='Symbol',
            value=str(self.marker_symbol)
        )
        self.widg_updatecolor_button = widgets.Button(
            description='Update colors',
            layout=widgets.Layout(width='150px')
        )
        self.widg_reset_button = widgets.Button(
            description='Reset symbols',
            layout=widgets.Layout(width='150px')
        )
        self.widg_plot_name = widgets.Text(
            placeholder='plot',
            value='plot',
            description='Name',
            layout=widgets.Layout(width='150px')
        )
        self.widg_plot_format = widgets.Text(
            placeholder='png',
            value='png',
            description='Format',
            layout=widgets.Layout(width='150px')
        )
        self.widg_scale = widgets.Text(
            placeholder='1',
            value='1',
            description="Scale",
            layout=widgets.Layout(width='150px')
        )
        self.widg_print_button = widgets.Button(
            description='Print',
            layout=widgets.Layout(width='100px')
        )
        self.widg_print_out = widgets.Output(
            layout=widgets.Layout(width='300px')
        )
        self.widg_description = widgets.Label(
            value='Tick the box next to the cross symbols in order to choose which windows visualizes the next '
                  'structure selected in the map above.'
        )
        self.widg_colordescription = widgets.Label(
            value='In the boxes below, the colors used in the plot. Colors can be written as a text string, i.e. red, '
                  'green,...,  or in a rgb/a, hex format. '
        )
        self.widg_colordescription2 = widgets.Label(
            value="After modifying a specific field, click on the 'Update colors' button to display the changes in "
                  "the plot."
        )
        self.widg_printdescription = widgets.Label(
            value="Click 'Print' to export the plot in the desired format. The resolution of the image can be increased"
                  " by increasing the 'Scale' value."
        )
        file1 = open("./assets/compressed_sensing/cross.png", "rb")
        image1 = file1.read()
        self.widg_img1 = widgets.Image(
            value=image1,
            format='png',
            width=30,
            height=30,
        )
        file2 = open("./assets/compressed_sensing/cross2.png", "rb")
        image2 = file2.read()
        self.widg_img2 = widgets.Image(
            value=image2,
            format='png',
            width=30,
            height=30,
        )

    def f_x(self, x):
        # Gives the classifications line
        if self.current_features[0] == self.current_features[1]:
            return x
        else:
            return -x * self.coefficients[self.current_features[0]] / self.coefficients[self.current_features[1]] - \
                   self.intercept / self.coefficients[self.current_features[1]]

    def update_markers(self):
        # Markers size and symbol are updated simultaneously
        with self.fig.batch_update():
            self.scatter_RS.marker.size = self.RS_sizes
            self.scatter_ZB.marker.size = self.ZB_sizes
            self.scatter_RS.marker.symbol = self.RS_symbols
            self.scatter_ZB.marker.symbol = self.ZB_symbols

    def set_markers_size(self, feature='Default size'):
        # Defines the size of the markers based on the input feature.
        # In case of default feature all markers have the same size.
        # Points marked with x/cross are set with a specific size

        if feature == 'Default size':

            sizes_RS = [self.marker_size] * self.RS_npoints
            sizes_ZB = [self.marker_size] * self.ZB_npoints
            symbols_RS = self.RS_symbols
            symbols_ZB = self.ZB_symbols

            try:
                point = symbols_RS.index('x')
                sizes_RS[point] = self.cross_size
            except:
                try:
                    point = symbols_ZB.index('x')
                    sizes_ZB[point] = self.cross_size
                except:
                    pass
            try:
                point = symbols_RS.index('cross')
                sizes_RS[point] = self.cross_size
            except:
                try:
                    point = symbols_ZB.index('cross')
                    sizes_ZB[point] = self.cross_size
                except:
                    pass

            self.RS_sizes = sizes_RS
            self.ZB_sizes = sizes_ZB
        else:

            min_value = min(min(self.D_selected_df.loc[self.D_selected_df['Structure'] == 'RS'][feature]),
                            min(self.D_selected_df.loc[self.D_selected_df['Structure'] == 'ZB'][feature]))
            max_value = max(max(self.D_selected_df.loc[self.D_selected_df['Structure'] == 'RS'][feature]),
                            max(self.D_selected_df.loc[self.D_selected_df['Structure'] == 'ZB'][feature]))
            coeff = 2 * self.marker_size / (max_value - min_value)
            sizes_RS = self.marker_size / 2 + coeff * self.D_selected_df.loc[self.D_selected_df['Structure'] == 'RS'][
                feature]
            sizes_ZB = self.marker_size / 2 + coeff * self.D_selected_df.loc[self.D_selected_df['Structure'] == 'ZB'][
                feature]
            self.RS_sizes = sizes_RS
            self.ZB_sizes = sizes_ZB

    def handle_xfeat_change(self, change):
        # changes the feature plotted on the x-axis
        # separating line is modified accordingly
        self.fig.update_layout(
            xaxis_title=change.new,
        )
        self.current_features[0] = self.features.index(change.new)
        self.scatter_RS['x'] = self.D_selected_df.loc[self.D_selected_df['Structure'] == 'RS'][change.new].to_numpy()
        self.scatter_ZB['x'] = self.D_selected_df.loc[self.D_selected_df['Structure'] == 'ZB'][change.new].to_numpy()
        line_x = np.linspace(self.D_selected_df[change.new].min(), self.D_selected_df[change.new].max(), 1000)
        line_y = self.f_x(line_x)
        self.scatter_line['x'] = line_x
        self.scatter_line['y'] = line_y
        min_x = min(min(self.scatter_RS['x']), min(self.scatter_ZB['x']))
        max_x = max(max(self.scatter_RS['x']), max(self.scatter_ZB['x']))
        min_delta = 0.05 * abs(max_x - min_x)
        self.fig.layout['xaxis'].range = [min_x - min_delta, max_x + min_delta]

    def handle_yfeat_change(self, change):
        # changes the feature plotted on the x-axis
        # separating line is modified accordingly
        self.fig.update_layout(
            yaxis_title=change.new,
        )
        self.current_features[1] = self.features.index(change.new)
        self.scatter_RS['y'] = self.D_selected_df.loc[self.D_selected_df['Structure'] == 'RS'][change.new].to_numpy()
        self.scatter_ZB['y'] = self.D_selected_df.loc[self.D_selected_df['Structure'] == 'ZB'][change.new].to_numpy()
        line_x = np.linspace(self.D_selected_df[self.features[self.current_features[0]]].min(),
                             self.D_selected_df[self.features[self.current_features[0]]].max(), 1000)
        line_y = self.f_x(self.line_x)
        self.scatter_line['x'] = line_x
        self.scatter_line['y'] = line_y
        min_y = min(min(self.scatter_RS['y']), min(self.scatter_ZB['y']))
        max_y = max(max(self.scatter_RS['y']), max(self.scatter_ZB['y']))
        min_delta = 0.05 * abs(max_y - min_y)
        self.fig.layout['yaxis'].range = [min_y - min_delta, max_y + min_delta]

    def handle_markerfeat_change(self, change):
        self.set_markers_size(feature=change.new)
        self.update_markers()

    def display_button_l_clicked(self, button):

        # Actions are performed only if the string inserted in the text widget corresponds to an existing compound
        if self.widg_compound_text_l.value in self.D_selected_df['Chem Formula'].tolist():
            structure_l = self.D_selected_df[self.D_selected_df['Chem Formula'] ==
                                        self.widg_compound_text_l.value]['Structure'].values[0]
            self.viewer_l.script(
                "load data/compressed_sensing/structures/" + structure_l + "_structures/"
                + self.widg_compound_text_l.value + ".xyz")

            symbols_RS = self.RS_symbols
            symbols_ZB = self.ZB_symbols

            try:
                point = symbols_RS.index('x')
                symbols_RS[point] = self.marker_symbol
            except:
                try:
                    point = symbols_ZB.index('x')
                    symbols_ZB[point] = self.marker_symbol
                except:
                    pass
            if structure_l == 'RS':
                point = np.where(self.scatter_RS['text'] == self.widg_compound_text_l.value)[0][0]
                symbols_RS[point] = 'x'
            if structure_l == 'ZB':
                point = np.where(self.scatter_ZB['text'] == self.widg_compound_text_l.value)[0][0]
                symbols_ZB[point] = 'x'

            self.ZB_symbols = symbols_ZB
            self.RS_symbols = symbols_RS
            self.set_markers_size(feature=self.widg_featmarker.value)
            self.update_markers()

    def display_button_r_clicked(self, button):

        # Actions are performed only if the string inserted in the text widget corresponds to an existing compound
        if self.widg_compound_text_r.value in self.D_selected_df['Chem Formula'].tolist():
            structure_r = self.D_selected_df[self.D_selected_df['Chem Formula'] ==
                                        self.widg_compound_text_r.value]['Structure'].values[0]
            self.viewer_r.script(
                "load data/compressed_sensing/structures/" + structure_r + "_structures/"
                + self.widg_compound_text_r.value + ".xyz")

            symbols_RS = self.RS_symbols
            symbols_ZB = self.ZB_symbols

            try:
                point = symbols_RS.index('cross')
                symbols_RS[point] = self.marker_symbol
            except:
                try:
                    point = symbols_ZB.index('cross')
                    symbols_ZB[point] = self.marker_symbol
                except:
                    pass
            if structure_r == 'RS':
                point = np.where(self.scatter_RS['text'] == self.widg_compound_text_r.value)[0][0]
                symbols_RS[point] = 'cross'
            if structure_r == 'ZB':
                point = np.where(self.scatter_ZB['text'] == self.widg_compound_text_r.value)[0][0]
                symbols_ZB[point] = 'cross'

            self.RS_symbols = symbols_RS
            self.ZB_symbols = symbols_ZB
            self.set_markers_size(feature=self.widg_featmarker.value)
            self.update_markers()

    def updatecolor_button_clicked(self, button):

        try:
          self.scatter_RS.update(marker=dict(color=self.widg_rscolor.value))
        except:
           pass
        try:
          self.scatter_ZB.update(marker=dict(color=self.widg_zbcolor.value))
        except:
           pass
        try:
            self.fig.update_layout(plot_bgcolor=self.widg_bgcolor.value)
        except:
            pass

    def handle_fontfamily_change(self, change):

        self.fig.update_layout(
            font=dict(family=change.new)
        )

    def handle_fontsize_change(self, change):

        self.fig.update_layout(
            font=dict(size=change.new)
        )

    def handle_markersize_change(self, change):

        self.marker_size = int(change.new)
        self.set_markers_size(feature=self.widg_featmarker.value)
        self.update_markers()

    def handle_crossize_change(self, change):

        self.cross_size = int(change.new)
        self.set_markers_size(feature=self.widg_featmarker.value)
        self.update_markers()

    def handle_linewidth_change(self, change):

        self.line_width = change.new
        self.scatter_line.update(line=dict(width=int(self.widg_linewidth.value))),

    def handle_linestyle_change(self, change):

        self.scatter_line.update(line=dict(dash=change.new))

    def bgtoggle_button_clicked(self, button):

        if self.bg_toggle:
            self.bg_toggle = False
            self.fig.update_layout(
                plot_bgcolor='white'
            )
        else:
            self.bg_toggle = True
            self.fig.update_layout(
                plot_bgcolor=self.widg_bgcolor.value
            )

    def print_button_clicked(self, button):

        try:
            os.mkdir("./plots")
        except:
            pass
        file_name = self.widg_plot_name.value + '.' + self.widg_plot_format.value
        self.fig.write_image("./plots/" + file_name, scale=self.widg_scale.value)
        self.widg_print_out.clear_output()
        with self.widg_print_out:
            local_file = FileLink('./plots/'+file_name, result_html_prefix="Click here to download: ")
            display(local_file)
        # javascript("window.open('./plots/" + str(file_name) + "' )")

    def reset_button_clicked(self, button):

        self.RS_symbols = [self.marker_symbol] * self.RS_npoints
        self.ZB_symbols = [self.marker_symbol] * self.ZB_npoints
        self.set_markers_size(self.widg_featmarker.value)
        self.update_markers()

    def handle_checkbox_l(self, change):
        if change.new:
            self.widg_checkbox_r.value = False
        else:
            self.widg_checkbox_r.value = True

    def handle_checkbox_r(self, change):
        if change.new:
            self.widg_checkbox_l.value = False
        else:
            self.widg_checkbox_l.value = True

    def view_structure_RS_l(self, formula):
        self.viewer_l.script("load data/compressed_sensing/structures/RS_structures/" + formula + ".xyz")

    def view_structure_RS_r(self, formula):
        self.viewer_r.script("load data/compressed_sensing/structures/RS_structures/" + formula + ".xyz")

    def view_structure_ZB_l(self, formula):
        self.viewer_l.script("load data/compressed_sensing/structures/ZB_structures/" + formula + ".xyz")

    def view_structure_ZB_r(self, formula):
        self.viewer_r.script("load data/compressed_sensing/structures/ZB_structures/" + formula + ".xyz")

    def update_point_RS(self, trace, points, selector):
        # changes the points labeled with a cross on the map.
        if not points.point_inds:
            return

        symbols_RS = self.RS_symbols
        symbols_ZB = self.ZB_symbols

        # The element previously marked with x/cross is marked with circle as default value
        if self.widg_checkbox_l.value:
            try:
                point = symbols_RS.index('x')
                symbols_RS[point] = self.marker_symbol
            except:
                try:
                    point = symbols_ZB.index('x')
                    symbols_ZB[point] = self.marker_symbol
                except:
                    pass
        if self.widg_checkbox_r.value:
            try:
                point = symbols_RS.index('cross')
                symbols_RS[point] = self.marker_symbol
            except:
                try:
                    point = symbols_ZB.index('cross')
                    symbols_ZB[point] = self.marker_symbol
                except:
                    pass

        if self.widg_checkbox_l.value:
            symbols_RS[points.point_inds[0]] = 'x'
        if self.widg_checkbox_r.value:
            symbols_RS[points.point_inds[0]] = 'cross'

        self.RS_symbols = symbols_RS
        self.ZB_symbols = symbols_ZB
        self.set_markers_size(feature=self.widg_featmarker.value)
        self.update_markers()

        formula = trace['text'][points.point_inds[0]][0]
        if self.widg_checkbox_l.value:
            self.widg_compound_text_l.value = formula
            self.view_structure_RS_l(formula)
        if self.widg_checkbox_r.value:
            self.widg_compound_text_r.value = formula
            self.view_structure_RS_r(formula)

    def update_point_ZB(self, trace, points, selector):
        if not points.point_inds:
            return

        symbols_RS = self.RS_symbols
        symbols_ZB = self.ZB_symbols

        # The element previously marked with x/cross is marked with circle as default value
        if self.widg_checkbox_l.value:
            try:
                point = symbols_RS.index('x')
                symbols_RS[point] = self.marker_symbol
            except:
                try:
                    point = symbols_ZB.index('x')
                    symbols_ZB[point] = self.marker_symbol
                except:
                    pass
        if self.widg_checkbox_r.value:
            try:
                point = symbols_RS.index('cross')
                symbols_RS[point] = self.marker_symbol
            except:
                try:
                    point = symbols_ZB.index('cross')
                    symbols_ZB[point] = self.marker_symbol
                except:
                    pass

        if self.widg_checkbox_l.value:
            symbols_ZB[points.point_inds[0]] = 'x'
        if self.widg_checkbox_r.value:
            symbols_ZB[points.point_inds[0]] = 'cross'

        self.RS_symbols = symbols_RS
        self.ZB_symbols = symbols_ZB
        self.set_markers_size(feature=self.widg_featmarker.value)
        self.update_markers()

        formula = trace['text'][points.point_inds[0]][0]
        if self.widg_checkbox_l.value:
            self.widg_compound_text_l.value = formula
            self.view_structure_ZB_l(formula)
        if self.widg_checkbox_r.value:
            self.widg_compound_text_r.value = formula
            self.view_structure_ZB_r(formula)

    def show(self):

        self.widg_featx.observe(self.handle_xfeat_change, names='value')
        self.widg_featy.observe(self.handle_yfeat_change, names='value')
        self.widg_featmarker.observe(self.handle_markerfeat_change, names='value')
        self.widg_checkbox_l.observe(self.handle_checkbox_l, names='value')
        self.widg_checkbox_r.observe(self.handle_checkbox_r, names='value')
        self.widg_display_button_l.on_click(self.display_button_l_clicked)
        self.widg_display_button_r.on_click(self.display_button_r_clicked)
        self.widg_updatecolor_button.on_click(self.updatecolor_button_clicked)
        self.widg_reset_button.on_click(self.reset_button_clicked)
        self.widg_print_button.on_click(self.print_button_clicked)
        self.widg_bgtoggle_button.on_click(self.bgtoggle_button_clicked)
        self.widg_linestyle.observe(self.handle_linestyle_change, names='value')
        self.scatter_RS.on_click(self.update_point_RS)
        self.scatter_ZB.on_click(self.update_point_ZB)
        self.widg_markersize.observe(self.handle_markersize_change, names='value')
        self.widg_crosssize.observe(self.handle_crossize_change, names='value')
        self.widg_linewidth.observe(self.handle_linewidth_change, names='value')
        self.widg_fontfamily.observe(self.handle_fontfamily_change, names='value')
        self.widg_fontsize.observe(self.handle_fontsize_change, names='value')

        output_l = widgets.Output()
        output_r = widgets.Output()
        output_l.layout = widgets.Layout(width="400px", height='350px')
        output_r.layout = widgets.Layout(width="400px", height='350px')

        with output_l:
            display(self.viewer_l)
        with output_r:
            display(self.viewer_r)

        box_print = widgets.HBox([self.widg_plot_name, self.widg_plot_format, self.widg_scale,
                                  self.widg_print_button, self.widg_print_out])

        box_features = widgets.HBox([self.widg_featx, self.widg_featy, self.widg_featmarker])
        container = widgets.VBox([self.widg_printdescription, box_print, box_features, self.fig,
                                  self.widg_description,
                                  widgets.HBox([
                                      widgets.VBox(
                                          [widgets.HBox([self.widg_compound_text_l, self.widg_display_button_l,
                                                         self.widg_img1, self.widg_checkbox_l]),
                                           output_l]),
                                      widgets.VBox(
                                          [widgets.HBox([self.widg_compound_text_r, self.widg_display_button_r,
                                                         self.widg_img2, self.widg_checkbox_r]),
                                           output_r]),
                                      ])
                                  ])

        display(container)

    def plot_appearance(self):

     box = widgets.VBox([widgets.HBox([self.widg_markersize, self.widg_crosssize]),
                    widgets.HBox([self.widg_linewidth, self.widg_linestyle]),
                    widgets.HBox([self.widg_fontsize, self.widg_fontfamily]),
                    self.widg_colordescription, self.widg_colordescription2,
                    widgets.HBox([self.widg_rscolor, self.widg_zbcolor]),
                    widgets.HBox([self.widg_bgtoggle_button, self.widg_bgcolor]),
                    widgets.HBox([self.widg_updatecolor_button, self.widg_reset_button])])

     display(box)


