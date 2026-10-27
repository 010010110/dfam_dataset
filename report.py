from collections import Counter, defaultdict
import sys


def analyze_sanitized_fasta(fasta_path):
    sf_counts = Counter()
    order_counts = Counter()
    class_counts = Counter()
    hierarchy = defaultdict(lambda: defaultdict(Counter))

    total_sequences = 0

    with open(fasta_path, "r") as infile:
        for line in infile:
            if not line.startswith(">"):
                continue

            total_sequences += 1
            header = line.strip()[1:]

            # Extrai os metadados taxonômicos após o '#'
            if "#" in header:
                tax_part = header.split("#", 1)[1].split()[0]
            else:
                tax_part = header.split()[0]

            nodes = tax_part.split("/")

            main_class = nodes[0] if len(nodes) > 0 else "Indefinido"
            order = nodes[1] if len(nodes) > 1 else "Indefinido"
            superfamily = nodes[2] if len(nodes) > 2 else f"[Apenas Ordem: {order}]"

            class_counts[main_class] += 1
            order_counts[order] += 1
            sf_counts[superfamily] += 1
            hierarchy[main_class][order][superfamily] += 1

    # --- Relatório no Terminal ---
    print("\n" + "=" * 65)
    print(f"RELATÓRIO DE COMPOSIÇÃO TAXONÔMICA: {fasta_path}")
    print("=" * 65)
    print(f"Total de sequências analisadas: {total_sequences}\n")

    print(f"{'Superfamília':<35} {'Qtd':>10} {'% do Total':>15}")
    print("-" * 65)
    for sf, count in sf_counts.most_common():
        pct = (count / total_sequences * 100) if total_sequences > 0 else 0
        print(f"{sf:<35} {count:>10} {pct:>14.2f}%")
    print("-" * 65)

    print("\n--- RESUMO POR NÍVEIS SUPERIORES ---")
    print("\n[Classes]")
    for c, count in class_counts.most_common():
        pct = (count / total_sequences * 100) if total_sequences > 0 else 0
        print(f"  - {c:<20} {count:>8} ({pct:.2f}%)")

    print("\n[Ordens]")
    for ord_name, count in order_counts.most_common():
        pct = (count / total_sequences * 100) if total_sequences > 0 else 0
        print(f"  - {ord_name:<20} {count:>8} ({pct:.2f}%)")

    print("\n--- DISTRIBUIÇÃO HIERÁRQUICA COMPLETA ---")
    for cl in sorted(hierarchy.keys()):
        print(f"\n{cl}")
        for ord_name in sorted(hierarchy[cl].keys()):
            ord_total = sum(hierarchy[cl][ord_name].values())
            print(f"  └── {ord_name} (Total: {ord_total})")
            for sf, count in hierarchy[cl][ord_name].most_common():
                print(f"        ├── {sf:<25}: {count}")
    print("\n" + "=" * 65 + "\n")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Uso: python3 report_taxa.py <arquivo_sanitizado.fasta>")
        sys.exit(1)

    analyze_sanitized_fasta(sys.argv[1])