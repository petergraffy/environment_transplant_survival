"""Render publication PDFs for visual inspection without changing the figures."""
from pathlib import Path
import json
import sys
import pypdfium2 as pdfium
from PIL import Image, ImageDraw, ImageStat

root = Path.cwd()
dest = root / "output/revision_2026_checklist/las_era_lung/visual_qa"
dest.mkdir(parents=True, exist_ok=True)
folders = [
    "output/figures/figure1_cox_baseline_timevarying",
    "output/figures/prior_year_and_timevarying_pollution_subgroup_cox/las_era_20261015",
    "output/figures/prior_year_pollution_quartile_aalen_johansen_simplified",
    "output/figures/cohort_flow_diagram",
    "output/figures/pollution_absolute_1yr_survival_maps",
    "output/figures/candidates_heart_lung_absolute_risk_abc",
    "output/figures/transplant_centers_by_organ_four_panel",
    "output/figures/waitlist_candidates_zcta_map",
    "output/figures/waitlist_geography_ab_maps",
    "output/figures/prior_year_initial_listing_cox_model_set",
    "output/figures/pollution_severity_associations",
    "output/revision_2026_checklist/aj_with_risk_tables",
    "output/revision_2026_checklist/baseline_diagnostics",
    "output/revision_2026_checklist/timevarying_diagnostics",
]
Image.MAX_IMAGE_PIXELS = None
records = []
thumbs = []
for folder in folders:
    for path in sorted((root / folder).glob("*.pdf")):
        if path.name in {"candidates_heart_lung_absolute_risk_abc_updated.pdf", "pollution_severity_association_forest.pdf"}:
            continue
        if len(sys.argv) > 1 and sys.argv[1] not in path.name:
            continue
        doc = pdfium.PdfDocument(str(path))
        assert len(doc) == 1, path
        page = doc[0]
        scale = min(1.5, 1400 / page.get_width(), 2000 / page.get_height())
        bitmap = page.render(scale=scale)
        im = bitmap.to_pil().convert("RGB")
        assert max(ImageStat.Stat(im).stddev) > 2, path
        name = f"{path.parent.name}__{path.stem}.png"
        im.save(dest / name)
        textpage = page.get_textpage()
        text = textpage.get_text_range()
        assert len(text.strip()) > 20, path
        textpage.close()
        png = path.with_suffix(".png")
        dimensions = None
        if png.exists():
            with Image.open(png) as original:
                dimensions = original.size
                assert min(dimensions) >= 1000, png
        records.append({"path": path.relative_to(root).as_posix(), "pages": 1,
                        "png_dimensions": dimensions, "nonblank": True,
                        "preview": (dest / name).relative_to(root).as_posix()})
        im.thumbnail((430, 620))
        thumbs.append((path.stem, im.copy()))
        bitmap.close()
        page.close()
        doc.close()
        print(path.name, flush=True)

for offset in range(0, len(thumbs), 12):
    batch = thumbs[offset:offset+12]
    sheet = Image.new("RGB", (1800, ((len(batch)+3)//4)*660), "white")
    draw = ImageDraw.Draw(sheet)
    for k, (label, im) in enumerate(batch):
        x, y = (k % 4)*450, (k//4)*660
        draw.text((x+5,y+5),label[:64],fill="black")
        sheet.paste(im,(x+5,y+30))
    sheet.save(dest / f"contact_{offset//12+1}.png")
(dest / "checks.json").write_text(json.dumps(records,indent=2),encoding="utf-8")
print(f"PASS: {len(records)} nonblank single-page PDFs rendered.")
