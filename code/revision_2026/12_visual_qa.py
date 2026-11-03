"""Render final revision PDFs and verify their paired scientific figures."""
from pathlib import Path
import json
import pypdfium2 as pdfium
from pypdf import PdfReader
from PIL import Image, ImageDraw, ImageStat

root = Path("output/revision_2026_checklist")
qa = root / "visual_qa"
qa.mkdir(exist_ok=True)
records = []
groups = {}
for folder in ("aj_with_risk_tables", "baseline_diagnostics", "timevarying_diagnostics", "svi_validation"):
    images = []
    for path in sorted((root / folder).glob("*.pdf")):
        reader = PdfReader(path)
        assert len(reader.pages) == 1, path
        assert len(reader.pages[0].extract_text().strip()) > 20, path
        png = path.with_suffix(".png")
        with Image.open(png) as original:
            assert min(original.size) >= 1000, png
            assert max(ImageStat.Stat(original.convert("RGB")).stddev) > 2, png
            dimensions = original.size
        document = pdfium.PdfDocument(str(path))
        page = document[0]
        bitmap = page.render(scale=0.9)
        rendered = bitmap.to_pil().convert("RGB")
        out = qa / f"{folder}_{path.stem}.png"
        rendered.save(out)
        assert max(ImageStat.Stat(rendered).stddev) > 2, path
        images.append((path.stem, rendered.copy()))
        records.append({"file": str(path), "pages": 1, "png_dimensions": dimensions,
                        "rendered_pdf": str(out), "nonblank": True})
        bitmap.close()
        page.close()
        document.close()
    groups[folder] = images

for folder, images in groups.items():
    cols = 4 if len(images) > 2 else 2
    cell_w = 520
    cell_h = 980 if folder == "aj_with_risk_tables" else 350
    rows = (len(images) + cols - 1) // cols
    sheet = Image.new("RGB", (cols * cell_w, rows * cell_h), "white")
    draw = ImageDraw.Draw(sheet)
    for i, (label, im) in enumerate(images):
        x, y = (i % cols) * cell_w, (i // cols) * cell_h
        draw.text((x + 8, y + 6), label, fill="black")
        im.thumbnail((cell_w - 12, cell_h - 30))
        sheet.paste(im, (x + 6, y + 26))
    sheet.save(qa / f"contact_{folder}.png")

assert len(records) == 34, len(records)
(qa / "checks.json").write_text(json.dumps(records, indent=2), encoding="utf-8")
print(f"PASS: {len(records)} PDFs rendered; paired PNGs are high-resolution and nonblank.")
