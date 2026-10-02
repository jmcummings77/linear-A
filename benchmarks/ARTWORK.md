# Report artwork

The report uses adapted SVG outputs generated with
[Book of Shapes](https://bookofshapes.com/?sort=popular), created by Nikolaj
Sokolowski. The fixed compositions in `report-art.svg` belong to the report's
visual design. They are decorative, separate from its measured charts and live
matrix visualizations.

## Sources

Exports were obtained on 2026-10-02 using the site's Download SVG controls.
Unlisted controls retain the site's defaults at that time.

| Report motif | Source pattern | Selected settings |
| --- | --- | --- |
| Cube emblem | [Iso Cube Wireframe](https://bookofshapes.com/patterns/iso-cube-wireframe) | Size 5, gap 1 |
| Woven grid | [Woven Grid](https://bookofshapes.com/patterns/woven_grid) | Columns 8, rows 8, steps 4, displacement 4.5, seed 7 |
| Flow lines | [Flow Lines](https://bookofshapes.com/patterns/flow_lines) | Default export |
| Spiral | [Spiral Morph](https://bookofshapes.com/patterns/spiral_morph) | Default export |
| Interference mesh | [Interference Mesh](https://bookofshapes.com/patterns/interference-mesh) | Default export |
| Voxel sphere | [Iso Sphere](https://bookofshapes.com/patterns/iso-sphere) | Dimension 8, gap 1 |

Adaptations recolor the exported geometry for the report's light and dark themes,
adjust line weights for small motifs, round coordinates to three decimal places,
combine adjacent strokes, and reuse repeated cube cells. The original export hashes are recorded in
`artwork-sources.json`. The website's source code and node graphs are not included.

## Reuse terms

The [Book of Shapes licence](https://bookofshapes.com/license), reviewed on
2026-10-02, permits generated graphics in personal and commercial projects,
including modifications, without requiring credit. It excludes redistribution
as standalone pattern packs, clipart sets, template libraries, or generators.
It also excludes three named recreations of existing artworks; none are used here.

These graphics are incorporated into the benchmark report. The repository's code
licence does not replace the artwork's separate terms. Follow the linked licence
when reusing the graphics beyond this report.

## Embedding

`report.py` validates and embeds the fixed SVG sprite directly in every generated
HTML file. Only passive geometry and local fragment references are permitted;
scripts, event handlers, external resources, and metadata are rejected. The
artwork needs no network access, does not animate, and is hidden from assistive
technology where nearby text already supplies the meaning.
