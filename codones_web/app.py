from flask import Flask, render_template, jsonify, request
from Bio import Entrez, SeqIO
from collections import Counter
import json

app = Flask(__name__)
Entrez.email = "204800@gmail.com"

# IDs de diferentes cepas de E. coli en GenBank (3+ opciones)
GENOMES = {
    "K-12 MG1655 (Referencia)": "NC_000913.3",
    "K-12 W3110": "NC_007779.1",
    "O157:H7 EDL933": "NC_002695.1",
    "CFT073 (UPEC)": "NC_004431.1",
}

GENOME_CACHE = {}

def descargar_genoma(genome_id):
    if genome_id in GENOME_CACHE:
        return GENOME_CACHE[genome_id]
    
    try:
        handle = Entrez.efetch(
            db="nucleotide",
            id=genome_id,
            rettype="gbwithparts",
            retmode="text"
        )
        registro = SeqIO.read(handle, "genbank")
        GENOME_CACHE[genome_id] = registro
        return registro
    except Exception as e:
        print(f"Error descargando {genome_id}: {e}")
        return None

def contenido_gc(secuencia):
    g = secuencia.count("G")
    c = secuencia.count("C")
    return (g + c) / len(secuencia) * 100 if len(secuencia) > 0 else 0

def contar_codones_genoma(secuencia):
    codones = [str(secuencia[i:i+3]) for i in range(0, len(secuencia) - 2)]
    return Counter(codones)

def extraer_genes_detallados(registro):
    """Extrae información detallada de todos los genes"""
    genes = []
    
    for idx, feature in enumerate(registro.features):
        if feature.type == "CDS":
            q = feature.qualifiers
            try:
                location_start = int(feature.location.start)
                location_end = int(feature.location.end)
                seq_gen = feature.extract(registro.seq)
                longitud = len(seq_gen)
                
                genes.append({
                    "id": idx,
                    "gen": q.get("gene", [f"Gene_{idx}"])[0],
                    "producto": q.get("product", ["Unknown"])[0],
                    "id_proteina": q.get("protein_id", ["-"])[0],
                    "longitud": longitud,
                    "inicio": location_start,
                    "fin": location_end,
                    "secuencia": str(seq_gen),
                    "gc_gen": contenido_gc(str(seq_gen)),
                    "nota": q.get("note", ["-"])[0] if "note" in q else "-"
                })
            except:
                continue
    
    return sorted(genes, key=lambda x: x['longitud'], reverse=True)

def analizar_espacios_basura(registro, genes_list):
    """Analiza espacios sin genes entre secuencias codificantes"""
    secuencia = str(registro.seq)
    espacios = []
    
    genes_ordenados = sorted(genes_list, key=lambda x: x['inicio'])
    
    for i in range(len(genes_ordenados) - 1):
        fin_actual = genes_ordenados[i]['fin']
        inicio_proximo = genes_ordenados[i+1]['inicio']
        
        if inicio_proximo > fin_actual:
            tamaño_espacio = inicio_proximo - fin_actual
            seq_espacio = secuencia[fin_actual:inicio_proximo]
            espacios.append({
                "entre_genes": f"{genes_ordenados[i]['gen']} -> {genes_ordenados[i+1]['gen']}",
                "posicion": f"{fin_actual}-{inicio_proximo}",
                "tamaño": tamaño_espacio,
                "contenido_gc": contenido_gc(seq_espacio),
                "secuencia_snippet": seq_espacio[:100] if len(seq_espacio) > 0 else ""
            })
    
    return sorted(espacios, key=lambda x: x['tamaño'], reverse=True)

def compactacion_genica(registro):
    """Calcula el porcentaje de compactación de genes"""
    genes = extraer_genes_detallados(registro)
    longitud_total = len(registro.seq)
    
    longitud_codificante = sum(g['longitud'] for g in genes)
    porcentaje_codificante = (longitud_codificante / longitud_total) * 100 if longitud_total > 0 else 0
    
    return {
        "longitud_total": longitud_total,
        "longitud_codificante": longitud_codificante,
        "porcentaje_compactacion": porcentaje_codificante,
        "numero_genes": len(genes),
        "densidad": len(genes) / (longitud_total / 1_000_000) if longitud_total > 0 else 0
    }

@app.route("/")
def index():
    registro = descargar_genoma(GENOMES["K-12 MG1655 (Referencia)"])
    if not registro:
        return "Error descargando genoma", 500
        
    secuencia_genoma = str(registro.seq).upper()
    
    genes = extraer_genes_detallados(registro)
    compactacion = compactacion_genica(registro)
    espacios = analizar_espacios_basura(registro, genes)
    
    gen_mas_grande = genes[0] if genes else None
    gen_mas_pequeno = genes[-1] if genes else None
    
    conteo_codones = contar_codones_genoma(secuencia_genoma)
    gc = contenido_gc(secuencia_genoma)
    
    codones_stop = {
        "TAA": conteo_codones.get("TAA", 0),
        "TAG": conteo_codones.get("TAG", 0),
        "TGA": conteo_codones.get("TGA", 0)
    }
    
    total_stop = sum(codones_stop.values())
    proporciones_stop = {k: (v / total_stop) * 100 if total_stop > 0 else 0 for k, v in codones_stop.items()}
    
    top_genes = genes[:15]
    
    # Obtener lista de genomas disponibles
    genomas_disponibles = list(GENOMES.keys())
    
    return render_template(
        "index.html",
        genome_length=len(secuencia_genoma),
        gc=gc,
        num_genes=len(genes),
        gene_density=compactacion["densidad"],
        compactacion=compactacion,
        gen_mas_grande=gen_mas_grande,
        gen_mas_pequeno=gen_mas_pequeno,
        top_genes=top_genes,
        stop_codons=codones_stop,
        stop_props=proporciones_stop,
        espacios=espacios[:10],
        genomas_disponibles=genomas_disponibles
    )

@app.route("/api/genes")
def api_genes():
    """API para obtener todos los genes con búsqueda"""
    registro = descargar_genoma(GENOMES["K-12 MG1655 (Referencia)"])
    if not registro:
        return jsonify({"error": "Error loading genome"}), 500
        
    genes = extraer_genes_detallados(registro)
    
    search = request.args.get('search', '').lower()
    if search:
        genes = [g for g in genes if search in g['gen'].lower() or search in g['producto'].lower()]
    
    return jsonify({
        "total": len(genes),
        "genes": genes[:100]
    })

@app.route("/api/gen/<int:gen_id>")
def api_gen_detalle(gen_id):
    """API para ver detalles completos de un gen específico"""
    registro = descargar_genoma(GENOMES["K-12 MG1655 (Referencia)"])
    if not registro:
        return jsonify({"error": "Error loading genome"}), 500
        
    genes = extraer_genes_detallados(registro)
    
    if gen_id >= len(genes):
        return jsonify({"error": "Gene not found"}), 404
    
    gen = genes[gen_id]
    secuencia = gen['secuencia']
    
    codones_gen = [secuencia[i:i+3] for i in range(0, len(secuencia) - 2, 3)]
    contador_codones = Counter(codones_gen)
    
    tabla_genetica = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    proteina = ''.join([tabla_genetica.get(codones_gen[i], '?') for i in range(len(codones_gen))])
    
    return jsonify({
        "gen": gen,
        "num_codones": len(codones_gen),
        "codones_frecuentes": dict(contador_codones.most_common(10)),
        "proteina_sintesis": proteina,
        "tamaño_proteina": len(proteina),
        "codones_detalles": list(contador_codones.items())
    })

@app.route("/api/comparar-genomas")
def api_comparar_genomas():
    """API para comparar 2 genomas específicos"""
    genoma1_nombre = request.args.get('genoma1')
    genoma2_nombre = request.args.get('genoma2')
    
    if not genoma1_nombre or not genoma2_nombre:
        return jsonify({"error": "Especifica genoma1 y genoma2"}), 400
    
    if genoma1_nombre not in GENOMES or genoma2_nombre not in GENOMES:
        return jsonify({"error": "Genoma no encontrado"}), 404
    
    reg1 = descargar_genoma(GENOMES[genoma1_nombre])
    reg2 = descargar_genoma(GENOMES[genoma2_nombre])
    
    if not reg1 or not reg2:
        return jsonify({"error": "Error descargando genomas"}), 500
    
    seq1 = str(reg1.seq).upper()
    seq2 = str(reg2.seq).upper()
    
    codones1 = contar_codones_genoma(seq1)
    codones2 = contar_codones_genoma(seq2)
    
    genes1 = extraer_genes_detallados(reg1)
    genes2 = extraer_genes_detallados(reg2)
    
    # Encontrar diferencias en codones
    todos_codones = set(codones1.keys()) | set(codones2.keys())
    diferencias_codones = []
    
    for codon in sorted(todos_codones):
        count1 = codones1.get(codon, 0)
        count2 = codones2.get(codon, 0)
        
        if count1 != count2:
            diff = count2 - count1
            diferencias_codones.append({
                "codon": codon,
                f"{genoma1_nombre}": count1,
                f"{genoma2_nombre}": count2,
                "diferencia": diff,
                "porcentaje_diff": (diff / count1 * 100) if count1 > 0 else 0
            })
    
    # Ordenar por mayor diferencia
    diferencias_codones.sort(key=lambda x: abs(x['diferencia']), reverse=True)
    
    comparativa = {
        "genoma1": {
            "nombre": genoma1_nombre,
            "longitud": len(seq1),
            "num_genes": len(genes1),
            "gc": contenido_gc(seq1),
            "compactacion": (sum(g['longitud'] for g in genes1) / len(seq1)) * 100,
            "densidad_atg": codones1.get("ATG", 0) / (len(seq1) / 1000) if len(seq1) > 0 else 0,
            "codones_top_10": dict(codones1.most_common(10)),
            "genes_top_5": [{"gen": g['gen'], "longitud": g['longitud'], "producto": g['producto'][:50]} for g in genes1[:5]]
        },
        "genoma2": {
            "nombre": genoma2_nombre,
            "longitud": len(seq2),
            "num_genes": len(genes2),
            "gc": contenido_gc(seq2),
            "compactacion": (sum(g['longitud'] for g in genes2) / len(seq2)) * 100,
            "densidad_atg": codones2.get("ATG", 0) / (len(seq2) / 1000) if len(seq2) > 0 else 0,
            "codones_top_10": dict(codones2.most_common(10)),
            "genes_top_5": [{"gen": g['gen'], "longitud": g['longitud'], "producto": g['producto'][:50]} for g in genes2[:5]]
        },
        "diferencias_codones": diferencias_codones[:30],  # Top 30 diferencias
        "resumen_diferencias": {
            "total_codones": len(todos_codones),
            "codones_diferentes": len(diferencias_codones),
            "diferencia_longitud": abs(len(seq2) - len(seq1)),
            "diferencia_genes": abs(len(genes2) - len(genes1)),
            "diferencia_gc": abs(contenido_gc(seq2) - contenido_gc(seq1))
        }
    }
    
    return jsonify(comparativa)

@app.route("/api/genomas-disponibles")
def api_genomas_disponibles():
    """Retorna lista de genomas disponibles"""
    return jsonify({
        "genomas": list(GENOMES.keys()),
        "total": len(GENOMES)
    })

@app.route("/api/genoma-completo")
def api_genoma_completo():
    """API para comparación detallada de todos los genomas"""
    registros = {}
    for nombre, gid in GENOMES.items():
        reg = descargar_genoma(gid)
        if reg:
            registros[nombre] = reg
    
    if not registros:
        return jsonify({"error": "No genomes loaded"}), 500
    
    comparativa = {}
    for nombre, registro in registros.items():
        secuencia = str(registro.seq).upper()
        genes = extraer_genes_detallados(registro)
        conteo = contar_codones_genoma(secuencia)
        
        comparativa[nombre] = {
            "longitud_genoma": len(secuencia),
            "num_genes": len(genes),
            "contenido_gc": contenido_gc(secuencia),
            "compactacion": (sum(g['longitud'] for g in genes) / len(secuencia)) * 100 if len(secuencia) > 0 else 0,
            "codones_top_20": dict(conteo.most_common(20)),
            "genes_top_10": [{"gen": g['gen'], "longitud": g['longitud'], "producto": g['producto']} for g in genes[:10]],
            "densidad_atg": conteo.get("ATG", 0) / (len(secuencia) / 1000) if len(secuencia) > 0 else 0
        }
    
    return jsonify(comparativa)

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=False)
