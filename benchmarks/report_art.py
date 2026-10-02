"""Load the report's fixed artwork as passive, self-contained SVG geometry."""
from pathlib import Path
import re
import xml.etree.ElementTree as ET

SVG_NAMESPACE = "http://www.w3.org/2000/svg"
COMMON = {"id", "fill", "stroke", "stroke-width", "stroke-linecap", "stroke-linejoin",
          "vector-effect", "transform"}
ATTRIBUTES = {
    "svg": {"class", "width", "height", "aria-hidden", "focusable"},
    "defs": set(), "symbol": {"id", "viewBox"}, "g": COMMON,
    "path": COMMON | {"d"}, "line": COMMON | {"x1", "y1", "x2", "y2"},
    "polygon": COMMON | {"points"}, "polyline": COMMON | {"points"},
    "rect": COMMON | {"x", "y", "width", "height", "rx", "ry"},
    "circle": COMMON | {"cx", "cy", "r"}, "ellipse": COMMON | {"cx", "cy", "rx", "ry"},
    "use": COMMON | {"href", "x", "y", "width", "height"},
}
VALUES = {
    "class": {"report-art-sprite"}, "aria-hidden": {"true"}, "focusable": {"false"},
    "fill": {"none", "currentColor", "var(--art-face-1)", "var(--art-face-2)"},
    "stroke": {"none", "currentColor"},
    "stroke-linecap": {"butt", "round", "square"},
    "stroke-linejoin": {"miter", "round", "bevel"},
    "vector-effect": {"non-scaling-stroke"},
}
IDENTIFIER = r"report-art-[a-z0-9-]+"
NUMBERS = r"[0-9eE+.,\s-]+"


def load_report_art(path=None):
    """Reject active content or outside resources before embedding the fixed sprite.

    The art is part of the report, not user-supplied benchmark data. Restricting
    its vocabulary keeps future artwork replacements offline and non-executable.
    Parsing and serializing also drops comments and processing instructions.
    """
    path = Path(path) if path is not None else Path(__file__).with_name("report-art.svg")
    text = path.read_text(encoding="utf-8")
    if len(text) > 500_000 or re.search(r"<!\s*(?:DOCTYPE|ENTITY)", text, re.I):
        raise ValueError("Report artwork must be a bounded SVG without document type declarations.")
    try:
        root = ET.fromstring(text)
    except ET.ParseError as error:
        raise ValueError("Report artwork is not valid XML.") from error
    if root.tag != f"{{{SVG_NAMESPACE}}}svg" or root.attrib != {
            "class": "report-art-sprite", "width": "0", "height": "0",
            "aria-hidden": "true", "focusable": "false"}:
        raise ValueError("Report artwork must be a hidden, nonfocusable SVG sprite.")

    identifiers, references = set(), set()
    for element in root.iter():
        if not element.tag.startswith(f"{{{SVG_NAMESPACE}}}"):
            raise ValueError("Report artwork contains an unsupported namespace.")
        tag = element.tag.split("}", 1)[1]
        if tag not in ATTRIBUTES or tag == "svg" and element is not root:
            raise ValueError("Report artwork contains unsupported elements.")
        if (element.text or "").strip() or (element.tail or "").strip():
            raise ValueError("Report artwork must contain only geometry.")
        for name, value in element.attrib.items():
            if name not in ATTRIBUTES[tag]:
                raise ValueError("Report artwork contains unsupported attributes.")
            if name in VALUES:
                valid = value in VALUES[name]
            elif name == "id":
                valid = re.fullmatch(IDENTIFIER, value) and value not in identifiers
                identifiers.add(value)
            elif name == "href":
                valid = re.fullmatch("#" + IDENTIFIER, value)
                references.add(value[1:])
            elif name == "d":
                valid = re.fullmatch(r"[MmLlHhVvCcSsQqTtAaZz0-9eE+.,\s-]+", value)
            elif name == "transform":
                valid = re.fullmatch(r"(?:\s*(?:matrix|translate|scale|rotate|skewX|skewY)\(" + NUMBERS + r"\)\s*)+", value)
            else:
                valid = re.fullmatch(NUMBERS, value)
            if not valid:
                raise ValueError("Report artwork contains an invalid attribute value.")
    if not references <= identifiers:
        raise ValueError("Report artwork contains unresolved local references.")
    ET.register_namespace("", SVG_NAMESPACE)
    return ET.tostring(root, encoding="unicode")
