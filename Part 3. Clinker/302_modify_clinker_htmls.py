print("\n\n\n\t**Running Script 302**\n\n")


# NOTE: this script requires Playwright + Chromium (pip install playwright ; playwright install chromium).
# Chromium is opened automatically (headless, no window) to render each Clinker HTML, export its SVG,

#---------------------------------------------------------------------------
# 1. set up
#---------------------------------------------------------------------------
#import modules
import json
from pathlib import Path
import base64
import html
import re
from playwright.sync_api import sync_playwright

#define pre-existing folders
p4_folder = Path(__file__).resolve().parent
p4_outfolder = p4_folder / "output_folder"
XerC_clinker_folder = p4_outfolder / "HTMLs from XerC 1kb clusters"
XerD_clinker_folder = p4_outfolder / "HTMLs from XerD 1kb clusters"

#confirm these folders exist
folders_to_check = [p4_outfolder, XerC_clinker_folder, XerD_clinker_folder]
for folder_path in folders_to_check:
    if not folder_path.is_dir():
        print(f"\n\tfolder not found: {folder_path}\n")
        raise FileNotFoundError

#create new output folders
SVG_images_folder = p4_outfolder / "polished images of 3kb clinkers"
SVG_images_folder.mkdir(parents=True, exist_ok=True)

temporary_folder = p4_outfolder / "temporary files"
temporary_folder.mkdir(parents=True, exist_ok=True)

#define input files
XerC_clinker_html_files = sorted(XerC_clinker_folder.glob("*.html"))
print(f"\nnumber of XerC Clinker HTMLs found: {len(XerC_clinker_html_files)}\n")
XerD_clinker_html_files = sorted(XerD_clinker_folder.glob("*.html"))                      
print(f"\nnumber of XerD Clinker HTMLs found: {len(XerD_clinker_html_files)}\n")


#---------------------------------------------------------------------------
# 2. settings to define and apply later for the Clinker HTML - Claude wrote this for me after it helped me figure out how to edit the HTML and flip the REV plasmids, I'm just editing
#---------------------------------------------------------------------------


#FYI
#most of this section number 2 was written by Claude on my request, because I do not know how to modify HTMLs



clinker_settings = {
    "plot":    {"transitionDuration": 0,
                "scaleFactor": 50},
    "cluster": {"nameFontSize": 10,
                "lociFontSize": 0,
                "spacing": 10},
    "gene":    {"shape": {"bodyHeight": 9,
                          "tipHeight": 5,
                          "tipLength": 12}},}

expose_anchor = ",y={check:t=>y.checkDomain"
expose_replacement = ",_x=(window._h=h,window._d=d),y={check:t=>y.checkDomain"
plot_block_old = """
} else {
  plot(data)
}
"""
plot_block_new = """
} else {
  plot(data);
  data.clusters.forEach(cluster => {
    if (cluster.name.includes("_REV_") && !cluster.name.includes("_FOR_")) {
      cluster.loci.forEach(locus => window._h.flip(locus));
    }
  });
  window._d.updateConfig(__SETTINGS__);
  window._d.update();
}
""".replace("__SETTINGS__", json.dumps(clinker_settings))
extract_svg_js = """
() => {
  //remove the (invisible) hover boxes/handles so they don't pad the bounding box
  document.querySelectorAll("g.hover").forEach(e => e.remove());

  const svg = document.getElementById("root-svg");
  const g = svg.querySelector("g.clusterMapG");
  const b = g.getBBox();
  const pad = 10;

  const clone = svg.cloneNode(true);
  clone.querySelector("g.clusterMapG").removeAttribute("transform");   //drop the zoom/pan transform
  clone.removeAttribute("cursor");
  clone.setAttribute("viewBox", `${b.x - pad} ${b.y - pad} ${b.width + 2 * pad} ${b.height + 2 * pad}`);
  clone.setAttribute("width", Math.ceil(b.width + 2 * pad));
  clone.setAttribute("height", Math.ceil(b.height + 2 * pad));
  return new XMLSerializer().serializeToString(clone);
}
"""
#---------------------------------------------------------------------------
# 3. set up PDF summary outfiles
#---------------------------------------------------------------------------
summary_pdf_template = """
<html><head><style>
  @page { size: A4 landscape; margin: 10mm; }
  body { margin: 0; font-family: system-ui, -apple-system, "Helvetica Neue", Arial, sans-serif; }
  .page { height: 186mm; page-break-after: always; }
  .page:last-child { page-break-after: auto; }
  h1 { font-size: 12pt; margin: 0 0 3mm 0; height: 6mm; }
  img { width: 100%; height: 175mm; object-fit: contain; object-position: left top; }
</style></head><body>
__PAGES__
</body></html>
"""


def natural_key(path):
    """sort so that cluster_2 comes before cluster_10"""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", path.stem)]

#---------------------------------------------------------------------------
# 4. start to modify the htmls
#---------------------------------------------------------------------------
def process_html(html_file, label, page):
    """modify one clinker html (via a temp copy) and save its SVG; returns the SVG path, or None on failure"""


    svg_file = SVG_images_folder / f"{label}_{html_file.stem.lower()}.svg"     # e.g. XerC_cluster_1.svg
    temp_html = temporary_folder / f"{label}_{html_file.name}"

    text = html_file.read_text()          #original is only ever read, never modified

    if text.count(expose_anchor) != 1:
        print(f"\t!! {html_file.name}: expose anchor found {text.count(expose_anchor)} times (expected 1), skipping")
        return None
    text = text.replace(expose_anchor, expose_replacement, 1)

    if text.count(plot_block_old) != 1:
        print(f"\t!! {html_file.name}: plot block found {text.count(plot_block_old)} times (expected 1), skipping")
        return None
    text = text.replace(plot_block_old, plot_block_new, 1)

    temp_html.write_text(text)

    try:
        page_errors = []
        page.on("pageerror", lambda err: page_errors.append(str(err)))

        page.goto(temp_html.resolve().as_uri())
        page.wait_for_selector("svg#root-svg g.gene polygon", timeout=30000)
        page.wait_for_timeout(1000)       #let the final update() finish drawing

        if page_errors:
            print(f"\t!! {html_file.name}: javascript errors in page: {page_errors}")
            return None

        svg_text = page.evaluate(extract_svg_js)
        svg_file.write_text(svg_text, encoding="utf-8")
        print(f"\tsaved {svg_file.name}")
        return svg_file

    finally:
        temp_html.unlink(missing_ok=True)   #delete the temporary HTML


#---------------------------------------------------------------------------
# 6. finish making the PDFs
#---------------------------------------------------------------------------
def make_summary_pdf(label, svg_files, browser):
    """combine the SVGs of one group (XerC or XerD) into a single multi-page PDF in p4_outfolder"""

    if not svg_files:
        print(f"\n\tno {label} images were exported, so no {label} summary PDF was made")
        return

    pdf_file = p4_outfolder / f"{label}_all_clusters_summary.pdf"

    pages = []
    for svg_file in sorted(svg_files, key=natural_key):
        encoded = base64.b64encode(svg_file.read_bytes()).decode("ascii")
        pages.append(f'<div class="page"><h1>{html.escape(svg_file.stem)}</h1>'
                     f'<img src="data:image/svg+xml;base64,{encoded}"></div>')
    document = summary_pdf_template.replace("__PAGES__", "\n".join(pages))

    page = browser.new_page()
    try:
        page.set_content(document)
        page.wait_for_timeout(500)
        page.pdf(path=str(pdf_file), prefer_css_page_size=True, print_background=True)
        print(f"\n\tsaved {pdf_file.name} ({len(svg_files)} pages)")
    finally:
        page.close()



#---------------------------------------------------------------------------
# 6. finish exporting the svgs and pdfs
#---------------------------------------------------------------------------

failed = []
saved_svgs = {"XerC": [], "XerD": []}

try:
    with sync_playwright() as pw:
        browser = pw.chromium.launch(headless=True)

        for label, html_files in [("XerC", XerC_clinker_html_files),
                                  ("XerD", XerD_clinker_html_files)]:
            print(f"\nprocessing {label} clinkers...")
            for html_file in html_files:
                page = browser.new_page()          #fresh page per file
                try:
                    svg_file = process_html(html_file, label, page)
                except Exception as err:
                    print(f"\t!! {html_file.name}: {err}")
                    svg_file = None
                finally:
                    page.close()

                if svg_file is None:
                    failed.append(f"{label}/{html_file.name}")
                else:
                    saved_svgs[label].append(svg_file)
        print("\nmaking summary PDFs...")
        for label in ["XerC", "XerD"]:
            try:
                make_summary_pdf(label, saved_svgs[label], browser)
            except Exception as err:
                print(f"\t!! could not make {label} summary PDF: {err}")
                failed.append(f"{label} summary PDF")





        
        browser.close()

#---------------------------------------------------------------------------
# 5. Clean up and finish
#---------------------------------------------------------------------------
finally:
    #clean up the temporary folder if nothing is left in it
    for leftover in temporary_folder.glob("*"):
        leftover.unlink()
    temporary_folder.rmdir()


if failed:
    print(f"\n\t{len(failed)} file(s) failed:")
    for name in failed:
        print(f"\t\t{name}")
else:
    print("\n\tno file failures: no errors were encountered")

print("\n\n\t***Finished running script 302!!***\n\n")
