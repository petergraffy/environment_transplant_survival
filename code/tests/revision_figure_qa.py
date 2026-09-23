"""Create inspection thumbnails without modifying publication figures."""
from pathlib import Path
import json
from PIL import Image, ImageDraw, ImageStat
from pypdf import PdfReader

Image.MAX_IMAGE_PIXELS = None
run = Path("output/revision_runs/rolling_20260922")
qa = run / "figure_qa"
qa.mkdir(exist_ok=True)
threshold = (run / "snapshot_complete.txt").stat().st_mtime
paths = sorted(p for p in Path("output/figures").rglob("*.png") if p.stat().st_mtime > threshold)
records = []
tiles = []
for path in paths:
    with Image.open(path) as original:
        original.load()
        width, height = original.size
        thumb = original.convert("RGB")
        thumb.thumbnail((560, 470))
        variance = sum(ImageStat.Stat(thumb).var)
        assert width > 100 and height > 100 and variance > 1, str(path)
        tile = Image.new("RGB", (600, 550), "#dddddd")
        tile.paste(thumb, ((600 - thumb.width) // 2, 10))
        label = path.name
        label = "\n".join(label[i:i+75] for i in range(0, len(label), 75))
        ImageDraw.Draw(tile).text((10, 487), label, fill="black")
        tiles.append(tile)
        records.append({"path": str(path), "width": width, "height": height, "variance": variance})
for start in range(0, len(tiles), 9):
    subset = tiles[start:start+9]
    sheet = Image.new("RGB", (1800, 550 * ((len(subset)+2)//3)), "white")
    for i, tile in enumerate(subset):
        sheet.paste(tile, ((i % 3)*600, (i//3)*550))
    sheet.save(qa / f"contact_{start//9+1}.png")
(qa / "image_checks.json").write_text(json.dumps(records, indent=2))
print(f"Checked {len(records)} updated PNG figures; inspection sheets saved to {qa}")
pdf_records = []
for path in sorted(Path("output/figures").rglob("*.pdf")):
    if path.stat().st_mtime <= threshold:
        continue
    try:
        reader = PdfReader(path)
        assert len(reader.pages) > 0
        pdf_records.append({"path": str(path), "pages": len(reader.pages), "valid": True})
    except Exception as error:
        pdf_records.append({"path": str(path), "valid": False, "error": str(error)})
(qa / "pdf_checks.json").write_text(json.dumps(pdf_records, indent=2))
failures = [row for row in pdf_records if not row["valid"]]
print(f"Checked {len(pdf_records)} PDFs; failures: {len(failures)}")
for row in failures:
    print(row)
