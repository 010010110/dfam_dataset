import os
import re
import sys


def parse_header(header_line):
    line = header_line.strip()
    if not line.startswith(">"):
        return None, None

    content = line[1:].strip()

    # Separa ID e o restante
    if "#" in content:
        seq_id, meta = content.split("#", 1)
    else:
        parts = content.split(None, 1)
        seq_id = parts[0]
        meta = parts[1] if len(parts) > 1 else ""

    # Extrai o nome da subfamília/elemento (último token da linha)
    raw_tokens = meta.strip().split()
    subfamily_name = raw_tokens[-1] if raw_tokens else ""

    # Se meta tiver só o nome do elemento, a taxonomia está em branco
    tax_str = raw_tokens[0] if len(raw_tokens) > 1 else meta

    clean_tax = re.sub(r"[;_/\-]", " ", tax_str).upper()
    tokens = set(clean_tax.split())
    name_clean = re.sub(r"[;_/\-]", " ", subfamily_name).upper()
    name_tokens = set(name_clean.split())

    # Bloqueia ruídos fora do escopo do ClassifyTE (Crypton, SVA, snoRNA)
    if any(b in clean_tax for b in ["CRYPTON", "SVA", "SNORNA"]):
        return None, None

    main_class = "Unknown"
    order = "Unknown"
    superfamily = None

    # -------------------------------------------------------------
    # 1. RETROTRANSPOSONS (Classe I)
    # -------------------------------------------------------------

    # --- 1.1 LTR ---
    # Prioridade para Gypsy verdadeiro (mesmo se contiver Gypsy-ERV na árvore)
    if "GYPSY" in tokens and not any(
        v in clean_tax for v in ["RETROVIRIDAE", "ORTHORETROVIRINAE", "MALR"]
    ):
        main_class, order, superfamily = "ClassI", "LTR", "Gypsy"
    elif "COPIA" in tokens:
        main_class, order, superfamily = "ClassI", "LTR", "Copia"
    elif "BEL" in tokens or "PAO" in tokens:
        main_class, order, superfamily = "ClassI", "LTR", "Bel-Pao"

    # --- 1.2 DIRS / Ngaro ---
    elif "DIRS" in tokens or "NGARO" in tokens:
        main_class, order = "ClassI", "DIRS"
        superfamily = "Ngaro" if "NGARO" in tokens else "DIRS"

    # --- 1.5 SINE ---
    elif (
        "SINE" in tokens
        or any(t in tokens for t in ["ALU", "MIR", "B2", "B4"])
        or "SINE" in clean_tax
    ):
        main_class, order = "ClassI", "SINE"
        if (
            "7SL" in clean_tax
            or "7SL RNA" in clean_tax
            or any(t in tokens for t in ["ALU", "B1"])
        ):
            superfamily = "7SL"
        elif "5S" in clean_tax or "5S RNA" in clean_tax:
            superfamily = "5S"
        elif (
            "TRNA" in clean_tax
            or "TRNA PROMOTER" in clean_tax
            or any(t in tokens for t in ["MIR", "B2"])
        ):
            superfamily = "tRNA"
        else:
            superfamily = None

    # --- 1.4 LINE (Apenas LINEs autônomos) ---
    elif (
        "DEPENDENT" not in clean_tax
        and "DERIVED" not in clean_tax
        and "SINE" not in tokens
    ):
        if "L1" in tokens or "LINE1" in tokens:
            main_class, order, superfamily = "ClassI", "LINE", "L1"
        elif "R2" in tokens or "R2" in name_tokens:
            main_class, order, superfamily = "ClassI", "LINE", "R2"
        elif "RTE" in tokens or "RTE" in name_tokens:
            main_class, order, superfamily = "ClassI", "LINE", "RTE"
        elif "JOCKEY" in tokens or "JOCKEY" in name_tokens:
            main_class, order, superfamily = "ClassI", "LINE", "Jockey"
        elif "I ELEMENT" in clean_tax or "I" in name_tokens:
            main_class, order, superfamily = "ClassI", "LINE", "I"
        elif "LINE" in tokens:
            main_class, order, superfamily = "ClassI", "LINE", None

    # -------------------------------------------------------------
    # 2. DNA TRANSPOSONS (Classe II - TIRs)
    # -------------------------------------------------------------
    is_class_ii = "CLASS II" in clean_tax or "DNA TRANSPOSITION" in clean_tax
    if order == "Unknown" and is_class_ii:
        main_class, order = "ClassII", "DNA"

        if (
            any(
                t in tokens
                for t in ["TC1", "MARINER", "TCMAR", "TC2", "TIGGER"]
            )
            or "TC1 MARINER" in clean_tax
        ):
            superfamily = "Tc1-Mariner"
        elif "HAT" in tokens:
            superfamily = "hAT"
        elif any(t in tokens for t in ["MUTATOR", "MULE", "MUDR"]):
            superfamily = "Mutator"
        elif "MERLIN" in tokens or "MERLIN" in name_tokens:
            superfamily = "Merlin"
        elif "TRANSIB" in tokens or "TRANSIB" in name_tokens:
            superfamily = "Transib"
        elif "P ELEMENT" in clean_tax or "P" in name_tokens:
            superfamily = "P"
        elif (
            "PIGGYBAC" in tokens
            or "PIGGYBAC" in name_tokens
            or "PB" in name_tokens
        ):
            superfamily = "PiggyBac"
        elif any(t in tokens for t in ["HARBINGER", "PIF"]):
            superfamily = "PIF-Harbinger"
        elif any(t in tokens for t in ["CACTA", "ENSPM", "CMC"]):
            superfamily = "CACTA"
        elif "TIR" in tokens or "TRANSPOSASE" in tokens:
            superfamily = None

    # Se a ordem não foi identificada, descarta
    if order == "Unknown" or main_class == "Unknown":
        return None, None

    # Monta taxonomia padronizada
    if superfamily:
        tax_path = f"{main_class}/{order}/{superfamily}"
    else:
        tax_path = f"{main_class}/{order}"

    sanitized_header = f">{seq_id}#{tax_path} {subfamily_name}".strip()
    return sanitized_header, order


def sanitize_file(
    input_path,
    output_path,
    discarded_path="descartados_classifyte_headers.txt",
):
    kept_count = 0
    discarded_count = 0

    for p in [output_path, discarded_path]:
        d = os.path.dirname(p)
        if d:
            os.makedirs(d, exist_ok=True)

    with (
        open(input_path, "r") as infile,
        open(output_path, "w") as outfile,
        open(discarded_path, "w") as discfile,
    ):
        write_current = False
        for line in infile:
            if line.startswith(">"):
                sanitized_header, order = parse_header(line)

                if sanitized_header:
                    write_current = True
                    outfile.write(f"{sanitized_header}\n")
                    kept_count += 1
                else:
                    write_current = False
                    discarded_count += 1
                    discfile.write(line)
            elif write_current:
                outfile.write(line)

    print("\n--- RESUMO DA SANITIZAÇÃO ---")
    print(f"Sequências mantidas (FASTA final): {kept_count}")
    print(f"Sequências descartadas:             {discarded_count}")
    print(f"FASTA salvo em:                     '{output_path}'")
    print(f"Descartados salvos em:              '{discarded_path}'\n")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(
            "Uso: python3 sanitize_classifyte.py <entrada.fasta> <saida.fasta> [descartados.txt]"
        )
        sys.exit(1)

    disc_out = (
        sys.argv[3] if len(sys.argv) > 3 else "descartados_classifyte_headers.txt"
    )
    sanitize_file(sys.argv[1], sys.argv[2], disc_out)