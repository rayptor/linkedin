from bokeh.plotting import figure, curdoc
from bokeh.models import ColumnDataSource, Div, Button, TextInput, Span, Spacer
from bokeh.layouts import column, row

def del_comma(x: int) -> str:
    return f"{x:,}".replace(",", "")

def collatz_max_altitude(n: int) -> int:
    peak = n
    x = n
    while x != 1:
        x = x >> 1 if x & 1 == 0 else 3 * x + 1
        if x > peak:
            peak = x
    return peak

def search_max_altitude(a: int, b: int) -> tuple[int, int]:
    max_n, max_alt = a, 0
    for n in range(a, b + 1):
        alt = collatz_max_altitude(n)
        if alt > max_alt:
            max_n, max_alt = n, alt
    return max_n, max_alt

def collatz(n: int) -> list:
    seq = [n]
    while n != 1:
        n = n >> 1 if n % 2 == 0 else 3 * n + 1
        seq.append(n)
    return seq

curdoc().theme = "dark_minimal"

p = figure(
    height=500, width=1000,
    x_axis_label = "Steps",
    y_axis_label = "Value (log scale)",
    y_axis_type="log",
    tools = "pan,wheel_zoom,box_zoom,save, reset",
    toolbar_location="above",
)
source = ColumnDataSource(data=dict(x=[], y=[]))
source_marker = ColumnDataSource(data=dict(x=[], y=[]))
r_line = p.line("x", "y", source=source, line_width=1.5, alpha=1.0)
r_circles = p.scatter("x", "y", source=source, size=4, alpha=0.5)
r_peak = p.scatter("x", "y",
                   source=source_marker, size=12, color="#ff0000",
                   marker="star_dot", line_color="#ff0000",
                   legend_label="Maximum altitude")

p.add_layout(Span(location=1, dimension='width',
                  line_color='#aaaaaa',
                  line_dash='dashed',
                  line_width=1))
p.legend.location = "top_right"
p.legend.label_text_font_size = "11px"

title_div = Div(
    text = "<h2 style='margin:0;font-family:sans-serif;color:#20262b'>"
         "Collatz: maximum altitude on [a,b]</h2>", width=1000
)

input_a = TextInput(title="", value="1", width=200)
input_b = TextInput(title="", value="1000", width=200)
row_a = row(Div(text="[a", width=8, styles={"line-height": "35px"}), input_a)
row_b = row(Div(text="b]", width=8, styles={"line-height": "35px"}), input_b)
inputs = row(row_a, row_b)

button_search = Button(label="Search", width=150, height=30, button_type="success")

check_div = Div(
    text = "", width = 1000,
    styles = {"font-family": "monospace", "font-size": "13px"}
)

result_div = Div(
    text = "", width = 1000,
    styles = {
        "font-family": "monospace", "font-size": "13px",
        "background": "#d4ddff", "padding": "10px",
        "border-radius": "6px", "border": "2px solid #20262b",
        "min-height": "36px"
    }
)


def search() -> None:
    try:
        a = int(input_a.value.strip())
        b = int(input_b.value.strip())
        assert 1 <= a < b
    except (ValueError, AssertionError):
        check_div.text = "<span style='color:#ff0000'>A < B and must be greater than 0.</span>"
        return

    if b - a > 10**6:
        check_div.text = "<span style='color:#ff0000'>Interval too large (max 10^6).</span>"
        return

    max_n, max_alt = search_max_altitude(a, b)
    sequence = collatz(max_n)
    peak_step = sequence.index(max_alt)
    source.data = dict(x=list(range(len(sequence))), y=sequence)
    source_marker.data = dict(x=[peak_step], y=[max_alt])

    r_line.glyph.line_color = "mediumslateblue"
    r_circles.glyph.fill_color = "skyblue"
    r_circles.glyph.line_color = "skyblue"

    p.title.text = (
        f"Collatz Sequence, N = {del_comma(max_n)}, "
        f"({del_comma(len(sequence)-1)} steps, Max Altitude = {del_comma(max_alt)})"
    )

    result_div.text = (
        f"<b>Interval</b> : [{del_comma(a)}, {del_comma(b)}], &nbsp;"
        f"<b>Number</b> : {del_comma(max_n)}, &nbsp;"
        f"<b>Max Altitude</b> : {del_comma(max_alt)}, &nbsp;"
        f"<b>Steps</b> : {del_comma(len(sequence)-1)}, &nbsp;"
        f"<b>Steps (peak)</b> : {del_comma(peak_step)}"
    )

button_search.on_click(lambda: search())
controls = row(Spacer(width=50), inputs, Spacer(height=20), button_search, align="start")
curdoc().add_root(column(title_div, controls, Spacer(width=50), check_div, p, result_div,
                         sizing_mode="fixed"))
curdoc().title = "Collatz Max Altitude"
