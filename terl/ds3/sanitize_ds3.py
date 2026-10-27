import os
import re
import sys


def parse_header(header_line):
    line = header_line.strip()
    if not line.startswith(">"):
        return header_line

    content = line[1:].strip()

    # Separa ID e metadados
    if "#" in content:
        seq_id, meta = content.split("#", 1)
    else:
        parts = content.split(None, 1)
        seq_id = parts[0]
        meta = parts[1] if len(parts) > 1 else ""

    # Normaliza texto para busca insensível a maiúsculas/minúsculas
    clean_meta = re.sub(r"[;_/\-]", " ", meta).upper()
    tokens = set(clean_meta.split())

    # -------------------------------------------------------------
    # 1. CLASSE (Avaliar Class II primeiro para evitar falso positivo em Class I)
    # -------------------------------------------------------------
    main_class = "Unknown"
    if "CLASS II" in clean_meta or "DNA TRANSPOSITION" in clean_meta:
        main_class = "ClassII"
    elif (
        "CLASS I" in clean_meta
        or "RETROTRANSPOSITION" in clean_meta
        or "RETROPOSON" in clean_meta
    ):
        main_class = "ClassI"

    # -------------------------------------------------------------
    # 2. SUPERFAMÍLIA & ORDEM
    # -------------------------------------------------------------
    superfamily = None
    order = "Unknown"

    # --- Gypsy (avaliado com prioridade, excluindo linhagens virais de ERV) ---
    if "GYPSY" in tokens and not any(
        v in clean_meta for v in ["RETROVIRIDAE", "ORTHORETROVIRINAE", "MALR"]
    ):
        order = "LTR"
        superfamily = "Gypsy"

    # --- ERV (Retrovírus endógenos legítimos e elementos associados) ---
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
        order = "LTR"
        superfamily = "ERV"

    # --- Demais Retrotransposons LTR ---
    elif "COPIA" in tokens:
        order = "LTR"
        superfamily = "Copia"
    elif "BEL" in tokens or "PAO" in tokens or "BEL PAO" in clean_meta:
        order = "LTR"
        superfamily = "Bel-Pao"

    # --- SINE (trata superfamílias específicas e evita duplicar SINE/SINE) ---
    elif any(k in tokens for k in ["ALU", "B2", "B4", "MIR"]):
        order = "SINE"
        if "ALU" in tokens:
            superfamily = "Alu"
        elif "MIR" in tokens:
            superfamily = "MIR"
        elif "B2" in tokens:
            superfamily = "B2"
        elif "B4" in tokens:
            superfamily = "B4"
    elif "SINE" in tokens:
        order = "SINE"
        superfamily = None  # Para na ordem SINE sem duplicar
    
    # --- LINE (exclui quem for SINE ou apenas dependente de LINE) ---
    elif (
        ("L1" in tokens or "LINE1" in tokens)
        and "SINE" not in tokens
        and "DEPENDENT" not in clean_meta
    ) or ("LINE" in tokens and "DEPENDENT" not in clean_meta):
        order = "LINE"
        superfamily = "L1"

    # --- Transposons de DNA (Classe II) ---
    elif (
        any(k in tokens for k in ["TC1", "MARINER", "TCMAR", "TC2", "TIGGER"])
        or "TC1 MARINER" in clean_meta
    ):
        order = "DNA"
        superfamily = "Tc1-Mariner"
    elif "HAT" in tokens:
        order = "DNA"
        superfamily = "hAT"
    elif any(k in tokens for k in ["MUTATOR", "MULE"]):
        order = "DNA"
        superfamily = "Mutator"
    elif any(k in tokens for k in ["HARBINGER", "PIF"]):
        order = "DNA"
        superfamily = "PIF-Harbinger"
    elif any(k in tokens for k in ["CACTA", "ENSPM"]):
        order = "DNA"
        superfamily = "CACTA"

    # -------------------------------------------------------------
    # 3. MONTAGEM DO CABEÇALHO REPEATMASKER
    # -------------------------------------------------------------
    # Casos com superfamília específica definida
    if superfamily:
        tax_path = f"{main_class}/{order}/{superfamily}"
    # Casos válidos apenas no nível de ordem (ex: SINEs genéricos como LFSINE)
    elif order == "SINE":
        tax_path = f"{main_class}/{order}"
    # Elementos fora dos alvos caem em Unknown
    elif order != "Unknown":
        tax_path = f"{main_class}/{order}/Unknown"
    else:
        tax_path = f"{main_class}/Unknown/Unknown"

    raw_parts = meta.strip().split()
    subfamily_name = raw_parts[-1] if raw_parts else ""

    return f">{seq_id}#{tax_path} {subfamily_name}".strip()


def sanitize_file(
    input_path, output_path, discarded_path="./terl/ds3/descartados_headers.txt"
):
    kept_count = 0
    discarded_count = 0

    # Cria os diretórios de saída se não existirem
    for file_path in [output_path, discarded_path]:
        dir_name = os.path.dirname(file_path)
        if dir_name:
            os.makedirs(dir_name, exist_ok=True)

    with (
        open(input_path, "r") as infile,
        open(output_path, "w") as outfile,
        open(discarded_path, "w") as discfile,
    ):
        write_current = False
        for line in infile:
            if line.startswith(">"):
                sanitized_header = parse_header(line)

                # Descarta se houver qualquer nível indefinido
                if "Unknown" in sanitized_header:
                    write_current = False
                    discarded_count += 1
                    discfile.write(line)
                else:
                    write_current = True
                    outfile.write(f"{sanitized_header}\n")
                    kept_count += 1
            elif write_current:
                outfile.write(line)

    print("\n--- RESUMO DA SANITIZAÇÃO (DS3) ---")
    print(f"Sequências mantidas (FASTA final): {kept_count}")
    print(f"Sequências descartadas:             {discarded_count}")
    print(f"FASTA sanitizado salvo em:          '{output_path}'")
    print(f"Lista de descartados salva em:      '{discarded_path}'\n")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(
            "Uso: python3 sanitize_ds3.py <entrada.fasta> <saida.fasta> [descartados.txt]"
        )
        sys.exit(1)

    discarded_out = (
        sys.argv[3]
        if len(sys.argv) > 3
        else "./terl/ds3/descartados_headers.txt"
    )
    sanitize_file(sys.argv[1], sys.argv[2], discarded_out)