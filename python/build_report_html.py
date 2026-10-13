"""Render results/REPORT.md to results/report.html (figures referenced relatively)."""
import os
import re
import markdown

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")

CSS = """
<title>Cooling Cloud Perturbations</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Source+Serif+4:opsz,wght@8..60,400;8..60,600;8..60,700&family=Source+Sans+3:wght@400;600&family=JetBrains+Mono:wght@400&display=swap">
<style>
/* layout: single reading column, figures break out wider; numbered sections mirror the report */
:root {
  --bg: #f7f8fa; --panel: #ffffff; --fg: #15191f; --muted: #555d68;
  --rule: #d9dde3; --accent: #1f5fae; --code-bg: #eef1f5; --fig-bg: #fcfcfb;
  --serif: "Source Serif 4", Georgia, "Times New Roman", serif;
  --sans: "Source Sans 3", "Helvetica Neue", Arial, sans-serif;
  --mono: "JetBrains Mono", ui-monospace, Menlo, Consolas, monospace;
}
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg: #12151a; --panel: #1a1e25; --fg: #e6e9ee; --muted: #a3abb6;
  --rule: #2c323b; --accent: #7fb0f0; --code-bg: #222831; --fig-bg: #fcfcfb; color-scheme: dark } }
:root[data-theme="dark"] {
  --bg: #12151a; --panel: #1a1e25; --fg: #e6e9ee; --muted: #a3abb6;
  --rule: #2c323b; --accent: #7fb0f0; --code-bg: #222831; --fig-bg: #fcfcfb; color-scheme: dark }
body { background: var(--bg); color: var(--fg); font-family: var(--sans); font-size: 17px; line-height: 1.6; }
main { max-width: 780px; margin: 0 auto; padding-inline: 20px; padding-block: 40px 80px; }
h1, h2, h3 { font-family: var(--serif); text-wrap: balance; line-height: 1.2; }
h1 { font-size: 2.1rem; font-weight: 700; margin: 0 0 0.6rem; }
h2 { font-size: 1.45rem; font-weight: 600; margin: 2.6rem 0 0.8rem; padding-top: 1.2rem; border-top: 1px solid var(--rule); }
h3 { font-size: 1.15rem; font-weight: 600; margin: 1.8rem 0 0.5rem; }
p, li { max-width: 68ch; }
a { color: var(--accent); }
strong { font-weight: 600; }
code { font-family: var(--mono); font-size: 0.86em; background: var(--code-bg); padding: 0.1em 0.35em; border-radius: 3px; }
ol li, ul li { margin-bottom: 0.45rem; }
.tbl { overflow-x: auto; margin: 1rem 0; }
table { border-collapse: collapse; font-size: 0.92rem; font-variant-numeric: tabular-nums; }
th, td { border-bottom: 1px solid var(--rule); padding: 0.4rem 0.8rem; text-align: left; }
th { font-weight: 600; color: var(--muted); font-size: 0.82rem; letter-spacing: 0.02em; }
figure { margin: 1.4rem 0; }
figure img { display: block; width: 100%; max-width: 100%; height: auto; background: var(--fig-bg); border: 1px solid var(--rule); border-radius: 4px; }
@media (min-width: 1100px) { figure { margin-inline: -150px; } }
.lede { color: var(--muted); font-size: 1.02rem; }
.meta { font-family: var(--mono); font-size: 0.78rem; color: var(--muted); letter-spacing: 0.02em; }
</style>
"""


def main():
    md = open(os.path.join(RES, "REPORT.md")).read()
    html = markdown.markdown(md, extensions=["tables", "sane_lists"])
    # images -> figures
    html = re.sub(r'<p><img alt="([^"]*)" src="([^"]+)" ?/?></p>',
                  r'<figure><img alt="\1" src="\2" loading="lazy"></figure>', html)
    html = html.replace("<table>", '<div class="tbl"><table>').replace("</table>", "</table></div>")
    # first paragraph after h1 is the lede
    html = re.sub(r"(</h1>\s*)<p>", r'\1<p class="meta">1D gPLUTO survey · linear theory · Oct 2026</p><p class="lede">', html, count=1)
    out = CSS + "<main>\n" + html + "\n</main>\n"
    path = os.path.join(RES, "report.html")
    open(path, "w").write(out)
    print("wrote", path)


if __name__ == "__main__":
    main()
