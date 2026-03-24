#!/usr/bin/env python3
from __future__ import annotations
import re
import sys
from pathlib import Path

# ------------------------------------------------------------
# Uso:
#   python make_array_geo.py input.geo output.geo NX NY
# ------------------------------------------------------------

def die(msg: str) -> None:
    raise SystemExit(msg)

def extract_blocks(text: str):
    m_cell11 = re.search(r'// Cell \(1,1\)', text)
    m_cell21 = re.search(r'// Cell \(2,1\)', text)
    m_footer = re.search(r'\nCoherence;\n', text)

    if not m_cell11:
        die("No encontré '// Cell (1,1)' en el archivo.")
    if not m_cell21:
        die("No encontré '// Cell (2,1)' en el archivo.")
    if not m_footer:
        die("No encontré 'Coherence;' en el archivo.")

    header = text[:m_cell11.start()]
    cell_template = text[m_cell11.start():m_cell21.start()]
    footer = text[m_footer.start()+1:]  # incluye Coherence; y lo de abajo
    return header, cell_template, footer


def update_header(header: str, nx: int, ny: int) -> str:
    header = re.sub(r'NX\s*=\s*\d+\s*;', f'NX = {nx};', header)
    header = re.sub(r'NY\s*=\s*\d+\s*;', f'NY = {ny};', header)
    return header


# ------------------------------------------------------------
# IDs locales de una celda
# Ajusta estos grupos si tu bloque base tiene otras entidades.
# ------------------------------------------------------------
POINT_IDS = {
    1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
    101,102,103,104
}

CURVE_IDS = {
    1,2,3,4,5,6,7,8,
    21,22,23,24,25,26,27,28,
    31,32,33,34,35,36,37,38,
    41,42,43,44,45,46,47,48,
    51,52,53,54,55,56,57,58,
    72,74,76,78
}

SURFACE_IDS = {
    101,102,103,104,105,106,107,108,
    201,202,203,204,205,206,207,208,
    211,212,213,214
}


def cell_offset(ix: int, iy: int, nx: int) -> int:
    # mismo esquema que venías usando:
    # serial = ix + nx*iy
    # offset = 10000 * serial
    return 10000 * (ix + nx * iy)


def replace_only_ids(body: str, allowed: set[int], offset: int) -> str:
    def repl(m: re.Match[str]) -> str:
        val = int(m.group(0))
        return str(offset + val) if val in allowed else m.group(0)
    return re.sub(r'\b\d+\b', repl, body)


def replace_signed_curve_ids(body: str, offset: int) -> str:
    def repl(m: re.Match[str]) -> str:
        sign = '-' if m.group(1) else ''
        val = int(m.group(2))
        newv = offset + val if val in CURVE_IDS else val
        return f"{sign}{newv}"
    return re.sub(r'(-?)(\d+)', repl, body)


def transform_line(line: str, ix: int, iy: int, nx: int) -> str:
    offset = cell_offset(ix, iy, nx)
    out = line
    stripped = out.strip()

    # Comentario de celda
    out = out.replace("// Cell (1,1)", f"// Cell ({ix+1},{iy+1})")

    # Traslaciones de geometría del bloque base
    out = out.replace("(0)*dx", f"({ix})*dx")
    out = out.replace("(0)*dy", f"({iy})*dy")

    # ----------------------------
    # Point(id) = {...};
    # ----------------------------
    m = re.match(r'^Point\((\d+)\)\s*=\s*\{([^}]*)\};$', stripped)
    if m:
        eid = int(m.group(1))
        body = m.group(2)
        return f"Point({offset + eid}) = {{{body}}};"

    # ----------------------------
    # Line(id) = {p1,p2};
    # Circle(id) = {p1,p2,p3};  o variantes
    # ----------------------------
    m = re.match(r'^(Line|Circle)\((\d+)\)\s*=\s*\{([^}]*)\};$', stripped)
    if m:
        kind = m.group(1)
        eid = int(m.group(2))
        body = m.group(3)
        body_new = replace_only_ids(body, POINT_IDS, offset)
        return f"{kind}({offset + eid}) = {{{body_new}}};"

    # ----------------------------
    # Curve Loop(id) = {...};
    # ----------------------------
    m = re.match(r'^Curve Loop\((\d+)\)\s*=\s*\{([^}]*)\};$', stripped)
    if m:
        eid = int(m.group(1))
        body = m.group(2)
        body_new = replace_signed_curve_ids(body, offset)
        return f"Curve Loop({offset + eid}) = {{{body_new}}};"

    # ----------------------------
    # Plane Surface(id) = {...};
    # ----------------------------
    m = re.match(r'^Plane Surface\((\d+)\)\s*=\s*\{([^}]*)\};$', stripped)
    if m:
        eid = int(m.group(1))
        body = m.group(2)
        body_new = replace_only_ids(body, SURFACE_IDS, offset)
        return f"Plane Surface({offset + eid}) = {{{body_new}}};"

    # ----------------------------
    # Transfinite Line{...} = ...
    # ----------------------------
    if stripped.startswith("Transfinite Line"):
        def repl(m: re.Match[str]) -> str:
            ids = [x.strip() for x in m.group(1).split(",")]
            new_ids = []
            for x in ids:
                val = int(x)
                new_ids.append(str(offset + val))
            return "{" + ",".join(new_ids) + "}"
        return re.sub(r'\{([^}]*)\}', repl, out, count=1)

    # ----------------------------
    # Transfinite Surface{...};
    # Recombine Surface{...};
    # ----------------------------
    if stripped.startswith("Transfinite Surface") or stripped.startswith("Recombine Surface"):
        def repl(m: re.Match[str]) -> str:
            ids = [x.strip() for x in m.group(1).split(",")]
            new_ids = []
            for x in ids:
                val = int(x)
                new_ids.append(str(offset + val))
            return "{" + ",".join(new_ids) + "}"
        return re.sub(r'\{([^}]*)\}', repl, out, count=1)

    # ----------------------------
    # out[] = Extrude {0,0,th} { Surface{101}; Layers{nZ}; Recombine; };
    # ----------------------------
    m = re.match(
        r'^out\[\]\s*=\s*Extrude\s*\{0,0,th\}\s*\{\s*Surface\{(\d+)\};\s*Layers\{nZ\};\s*Recombine;\s*\};$',
        stripped
    )
    if m:
        sid = int(m.group(1))
        return f"out[] = Extrude {{0,0,th}} {{ Surface{{{offset + sid}}}; Layers{{nZ}}; Recombine; }};"

    # ----------------------------
    # slave_z[] += {236};
    # ----------------------------
    m = re.match(r'^slave_z\[\]\s*\+=\s*\{(\d+)\};$', stripped)
    if m:
        sid = int(m.group(1))
        return f"slave_z[] += {{{offset + sid}}};"

    return out


def make_cell(cell_template: str, ix: int, iy: int, nx: int) -> str:
    lines = cell_template.splitlines()
    new_lines = [transform_line(line, ix, iy, nx) for line in lines]
    return "\n".join(new_lines) + "\n"


def build_geo(input_geo: Path, output_geo: Path, nx: int, ny: int) -> None:
    text = input_geo.read_text()

    header, cell_template, footer = extract_blocks(text)
    header = update_header(header, nx, ny)

    cells = []
    for iy in range(ny):
        for ix in range(nx):
            cells.append(make_cell(cell_template, ix, iy, nx))

    new_text = header + "".join(cells) + footer
    output_geo.write_text(new_text)


def main():
    if len(sys.argv) != 5:
        die("Uso: python make_array_geo.py input.geo output.geo NX NY")

    input_geo = Path(sys.argv[1])
    output_geo = Path(sys.argv[2])
    nx = int(sys.argv[3])
    ny = int(sys.argv[4])

    if not input_geo.exists():
        die(f"No existe: {input_geo}")

    build_geo(input_geo, output_geo, nx, ny)
    print(f"Archivo generado: {output_geo}")


if __name__ == "__main__":
    main()