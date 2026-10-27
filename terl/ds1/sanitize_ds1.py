from collections import Counter
import os
import re
import sys


def classify_header(header_line):
    line = header_line.strip()
    if not line.startswith(">"):
        return None, None

    content = line[1:].strip()

    # Separa ID e metadados
    if "#" in content:
        seq_id, meta = content.split("#", 1)
    else:
        parts = content.split(None, 1)
        seq_id = parts[0]
        meta = parts[1] if len(parts) > 1 else ""

    # Normaliza texto separando termos por espaço
    clean_meta = re.sub(r"[;_/\-]", " ", meta).upper()
    tokens = set(clean_meta.split())

    resolved_class = None
    resolved_order = None
    resolved_superfamily = None

    # -------------------------------------------------------------
    # 1. Gypsy (avaliado com prioridade, excluindo quem for vírus/ERV)
    # -------------------------------------------------------------
    if "GYPSY" in tokens and not any(
        v in clean_meta for v in ["RETROVIRIDAE", "ORTHORETROVIRINAE", "MALR"]
    ):
        resolved_class, resolved_order, resolved_superfamily = (
            "ClassI",
            "LTR",
            "Gypsy",
        )

    # -------------------------------------------------------------
    # 2. ERV (Retrovírus endógenos legítimos)
    # -------------------------------------------------------------
    elif any(
        k in clean_meta
        for k in [
            "RETROVIRIDAE",
            "ORTHORETROVIRINAE",
            "MALR",
            "HERV",
            "ERV1",
            "ERV2",
            "ERV3",
            "ERVK",
            "ERVL",
        ]
    ) or ("ERV" in tokens and "GYPSY" not in tokens):
        resolved_class, resolved_order, resolved_superfamily = (
            "ClassI",
            "LTR",
            "ERV",
        )

    # -------------------------------------------------------------
    # 3. Copia
    # -------------------------------------------------------------
    elif "COPIA" in tokens:
        resolved_class, resolved_order, resolved_superfamily = (
            "ClassI",
            "LTR",
            "Copia",
        )

    # -------------------------------------------------------------
    # 4. Bel-Pao
    # -------------------------------------------------------------
    elif "BEL" in tokens or "PAO" in tokens or "BEL PAO" in clean_meta:
        resolved_class, resolved_order, resolved_superfamily = (
            "ClassI",
            "LTR",
            "Bel-Pao",
        )

    # -------------------------------------------------------------
    # 5. L1 (Genuíno - descarta pseudogenes/SINEs dependentes de L1)
    # -------------------------------------------------------------
    elif (
        (
            ("L1" in tokens or "LINE1" in tokens)
            and "DEPENDENT" not in clean_meta
            and "DERIVED" not in clean_meta
        )
        or ("LINE" in tokens and "DEPENDENT" not in clean_meta)
    ) and "SINE" not in tokens:
        resolved_class, resolved_order, resolved_superfamily = (
            "ClassI",
            "LINE",
            "L1",
        )

    # -------------------------------------------------------------
    # 6. Tc1-Mariner
    # -------------------------------------------------------------
    elif (
        any(k in tokens for k in ["TCMAR", "TC1", "MARINER", "TC2", "TIGGER"])
        or "TC1 MARINER" in clean_meta
    ):
        resolved_class, resolved_order, resolved_superfamily = (
            "ClassII",
            "DNA",
            "Tc1-Mariner",
        )

    # -------------------------------------------------------------
    # 7. hAT
    # -------------------------------------------------------------
    elif "HAT" in tokens:
        resolved_class, resolved_order, resolved_superfamily = (
            "ClassII",
            "DNA",
            "hAT",
        )

    # Se não casar estritamente com as 7 superfamílias, descarta
    if not resolved_superfamily:
        return None, None

    raw_tokens = meta.strip().split()
    subfamily_name = raw_tokens[-1] if raw_tokens else ""

    sanitized_header = f">{seq_id}#{resolved_class}/{resolved_order}/{resolved_superfamily} {subfamily_name}".strip()
    return sanitized_header, resolved_superfamily


def filter_and_sanitize(
    input_fasta, output_fasta, discarded_path="./terl/ds1/descartados_headers.txt"
):
    counts = Counter()
    total_in = 0
    total_kept = 0
    discarded_count = 0
    write_current = False

    for path in [output_fasta, discarded_path]:
        d = os.path.dirname(path)
        if d:
            os.makedirs(d, exist_ok=True)

    with (
        open(input_fasta, "r") as infile,
        open(output_fasta, "w") as outfile,
        open(discarded_path, "w") as discfile,
    ):
        for line in infile:
            if line.startswith(">"):
                total_in += 1
                new_header, sf = classify_header(line)

                if new_header:
                    write_current = True
                    outfile.write(f"{new_header}\n")
                    counts[sf] += 1
                    total_kept += 1
                else:
                    write_current = False
                    discarded_count += 1
                    discfile.write(line)
            elif write_current:
                outfile.write(line)

    print("\n--- RESUMO DO PROCESSAMENTO (DS1) ---")
    print(f"Total de sequências de entrada: {total_in}")
    print(f"Total mantidas (7 superfamílias): {total_kept}")
    print(f"Total descartadas:                {discarded_count}")
    print("\nContagem por Superfamília mantida:")
    for sf, qtd in sorted(counts.items()):
        print(f"  - {sf:15s}: {qtd} sequências")
    print(f"\nFASTA salvo em:                   '{output_fasta}'")
    print(f"Descartados salvos em:            '{discarded_path}'\n")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(
            "Uso: python3 sanitize_ds1.py <entrada.fasta> <saida.fasta> [descartados.txt]"
        )
        sys.exit(1)

    discarded_file = (
        sys.argv[3]
        if len(sys.argv) > 3
        else "./terl/ds1/descartados_headers.txt"
    )
    filter_and_sanitize(sys.argv[1], sys.argv[2], discarded_file)