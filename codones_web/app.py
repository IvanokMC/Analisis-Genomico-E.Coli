"""
Análisis Genómico de Escherichia coli K-12 MG1655
Proyecto de Bioinformática - Universidad Nacional San Antonio Abad del Cusco
"""

from flask import Flask, render_template, jsonify, request
from Bio import Entrez, SeqIO
from collections import Counter
import json
from groq import Groq

app = Flask(__name__)
Entrez.email = "204800@gmail.com"

# Configuración de Groq AI
GROQ_API_KEY = "inserte_api_key"
groq_client = Groq(api_key=GROQ_API_KEY)

# Contexto de conversación para el chatbot
CONVERSATION_CONTEXT = {}

# IDs de diferentes cepas de E. coli en GenBank
GENOMES = {
    "K-12 MG1655 (Referencia)": "NC_000913.3",
    "K-12 W3110": "NC_007779.1",
    "O157:H7 EDL933": "NC_002695.1",
    "CFT073 (UPEC)": "NC_004431.1",
    "BL21 (TaKaRa)": "CP010816.1",
}

# Caché en memoria para genomas descargados
GENOME_CACHE = {}

# Variable para la cepa seleccionada (por defecto K-12)
CEPA_ACTUAL = "NC_000913.3"

def descargar_genoma(genome_id):
    """
    Descarga un genoma desde NCBI GenBank.
    
    Args:
        genome_id (str): ID del genoma en GenBank (ej: "NC_000913.3")
    
    Returns:
        SeqRecord: Objeto BioPython con el genoma completo
    """
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
    """
    Calcula el porcentaje de contenido GC en una secuencia.
    
    Args:
        secuencia (str): Secuencia de ADN
    
    Returns:
        float: Porcentaje de G+C
    """
    g = secuencia.count("G")
    c = secuencia.count("C")
    return (g + c) / len(secuencia) * 100 if len(secuencia) > 0 else 0

def contar_tripletes_genoma_completo(secuencia):
    """
    Cuenta TODOS los tripletes en el genoma completo (análisis exploratorio).
    
    NOTA: Este método cuenta tripletes sin considerar marco de lectura ni
    si están dentro de genes. Es útil para análisis estadístico pero NO
    representa codones funcionales.
    
    Args:
        secuencia (str): Secuencia completa del genoma
    
    Returns:
        Counter: Diccionario con conteo de cada triplete
    """
    codones = [str(secuencia[i:i+3]) for i in range(0, len(secuencia) - 2)]
    return Counter(codones)

def contar_codones_en_cds(registro):
    """
    Cuenta codones SOLO dentro de CDS anotados (análisis biológicamente correcto).
    
    Este método:
    - Solo analiza secuencias codificantes (CDS)
    - Respeta el marco de lectura correcto (codones de 3 nucleótidos)
    - Cuenta codones funcionales reales
    
    Args:
        registro (SeqRecord): Registro de GenBank con anotaciones
    
    Returns:
        Counter: Diccionario con conteo de cada codón funcional
    """
    contador = Counter()
    
    for feature in registro.features:
        if feature.type == "CDS":
            secuencia_cds = str(feature.extract(registro.seq))
            # Dividir en codones respetando marco de lectura
            for i in range(0, len(secuencia_cds) - 2, 3):
                codon = secuencia_cds[i:i+3]
                if len(codon) == 3:
                    contador[codon] += 1
    
    return contador

def contar_stops_finales(registro):
    """
    Cuenta solo el último codón de cada gen (codones STOP reales).
    
    Este es el método más específico para validar que cada gen
    tiene exactamente un codón de terminación.
    
    Args:
        registro (SeqRecord): Registro de GenBank con anotaciones
    
    Returns:
        dict: {"TAA": count, "TAG": count, "TGA": count}
    """
    stops = {"TAA": 0, "TAG": 0, "TGA": 0}
    
    for feature in registro.features:
        if feature.type == "CDS":
            seq_cds = str(feature.extract(registro.seq))
            if len(seq_cds) >= 3:
                ultimo_codon = seq_cds[-3:]
                if ultimo_codon in stops:
                    stops[ultimo_codon] += 1
    
    return stops

def contar_atg_iniciales(registro):
    """
    Cuenta solo el primer codón (ATG) de cada gen (codones de inicio reales).
    
    Args:
        registro (SeqRecord): Registro de GenBank con anotaciones
    
    Returns:
        tuple: (contador_atg, lista_genes_sin_atg)
    """
    contador_atg = 0
    genes_sin_atg = []
    
    for feature in registro.features:
        if feature.type == "CDS":
            seq_cds = str(feature.extract(registro.seq))
            if len(seq_cds) >= 3:
                primer_codon = seq_cds[:3]
                if primer_codon == "ATG":
                    contador_atg += 1
                else:
                    # Registrar genes que no empiezan con ATG (casos raros)
                    gene_name = feature.qualifiers.get("gene", ["Unknown"])[0]
                    genes_sin_atg.append({
                        "gen": gene_name,
                        "primer_codon": primer_codon
                    })
    
    return contador_atg, genes_sin_atg

def extraer_genes_detallados(registro):
    """
    Extrae información detallada de todos los genes anotados.
    
    Args:
        registro (SeqRecord): Registro de GenBank
    
    Returns:
        list: Lista de diccionarios con información de cada gen
    """
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
                    "index_original": idx,  # NUEVO: Guardar índice original
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
    
    return genes  # NO ORDENAR AQUÍ

def analizar_espacios_basura(registro, genes_list):
    """
    Analiza espacios intergénicos (regiones no codificantes entre genes).
    
    Args:
        registro (SeqRecord): Registro de GenBank
        genes_list (list): Lista de genes extraídos
    
    Returns:
        list: Espacios intergénicos ordenados por tamaño
    """
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
    """
    Calcula estadísticas de compactación del genoma.
    
    Args:
        registro (SeqRecord): Registro de GenBank
    
    Returns:
        dict: Estadísticas de compactación génica
    """
    genes = extraer_genes_detallados(registro)
    longitud_total = len(registro.seq)
    
    longitud_codificante = sum(g['longitud'] for g in genes)
    porcentaje_codificante = (longitud_codificante / longitud_total) * 100 if longitud_total > 0 else 0
    
    return {
        "longitud_total": longitud_total,
        "longitud_codificante": longitud_codificante,
        "longitud_no_codificante": longitud_total - longitud_codificante,
        "porcentaje_compactacion": porcentaje_codificante,
        "numero_genes": len(genes),
        "densidad": len(genes) / (longitud_total / 1_000_000) if longitud_total > 0 else 0
    }

def validar_con_literatura(registro, stops_finales):
    """
    Valida resultados contra valores esperados de la literatura científica.
    
    Referencia: 
    - Blattner et al. (1997) "The complete genome sequence of Escherichia coli K-12"
    - Nakamura et al. (2000) "Codon usage tabulated from GenBank"
    
    Args:
        registro (SeqRecord): Registro de GenBank
        stops_finales (dict): Conteo de codones STOP finales
    
    Returns:
        dict: Resultados de validación
    """
    longitud = len(registro.seq)
    gc = contenido_gc(str(registro.seq).upper())
    
    genes = [f for f in registro.features if f.type == "CDS"]
    num_genes = len(genes)
    
    # Valores esperados de la literatura
    validaciones = {
        "longitud_ok": 4_500_000 <= longitud <= 4_700_000,
        "longitud_valor": longitud,
        "longitud_esperado": "4.5-4.7 Mb",
        
        "gc_ok": 50.5 <= gc <= 51.0,
        "gc_valor": round(gc, 2),
        "gc_esperado": "50.5-51.0%",
        
        "genes_ok": 4_200 <= num_genes <= 4_400,
        "genes_valor": num_genes,
        "genes_esperado": "4,200-4,400",
    }
    
    # Validar proporciones de codones STOP usando stops_finales
    codones_stop = stops_finales
    
    total_stop = sum(codones_stop.values())
    if total_stop > 0:
        prop_taa = (codones_stop["TAA"] / total_stop) * 100
        prop_tag = (codones_stop["TAG"] / total_stop) * 100
        prop_tga = (codones_stop["TGA"] / total_stop) * 100
        
        validaciones.update({
            "taa_ok": 60 <= prop_taa <= 65,
            "taa_valor": round(prop_taa, 2),
            "taa_esperado": "60-65%",
            
            "tag_ok": 5 <= prop_tag <= 10,
            "tag_valor": round(prop_tag, 2),
            "tag_esperado": "5-10%",
            
            "tga_ok": 25 <= prop_tga <= 35,
            "tga_valor": round(prop_tga, 2),
            "tga_esperado": "25-35%"
        })
    
    return validaciones

# ========================================
# SISTEMA DE IA - CHATBOT
# ========================================

def obtener_datos_genomicos_contexto():
    """Obtiene resumen de datos genómicos para contexto de IA."""
    try:
        registro = descargar_genoma(GENOMES["K-12 MG1655 (Referencia)"])
        if not registro:
            return "No se pudo cargar el genoma."
        
        secuencia = str(registro.seq).upper()
        genes = extraer_genes_detallados(registro)
        compactacion = compactacion_genica(registro)
        codones_cds = contar_codones_en_cds(registro)
        stops_finales = contar_stops_finales(registro)
        
        contexto = f"""
DATOS GENÓMICOS - E. coli K-12 MG1655:
- Longitud: {len(secuencia):,} pb ({len(secuencia)/1_000_000:.2f} Mb)
- Genes: {len(genes):,}
- GC: {contenido_gc(secuencia):.2f}%
- Densidad: {compactacion['porcentaje_compactacion']:.2f}%
- STOP: TAA={stops_finales['TAA']}, TAG={stops_finales['TAG']}, TGA={stops_finales['TGA']}
"""
        return contexto
    except Exception as e:
        return f"Error: {str(e)}"

def obtener_prompt_sistema(rol):
    """Genera el prompt del sistema según el rol."""
    contexto = obtener_datos_genomicos_contexto()
    
    config = {
        "estudiante": "Explicar de forma simple y didáctica. Usa analogías.",
        "investigador": "Análisis técnico profundo con terminología científica.",
        "docente": "Balance entre rigor científico y claridad pedagógica.",
        "divulgador": "Lenguaje accesible para audiencia general."
    }
    
    return f"""Eres un asistente de análisis genómico de E. coli.
ROL: {config.get(rol, config['estudiante'])}

{contexto}

Responde de forma clara y profesional."""

@app.route("/api/chat", methods=["POST"])
def api_chat():
    """Endpoint del chatbot IA."""
    try:
        data = request.get_json()
        mensaje = data.get("message", "")
        rol = data.get("role", "estudiante")
        session_id = data.get("session_id", "default")
        
        if not mensaje:
            return jsonify({"error": "Mensaje vacío"}), 400
        
        if session_id not in CONVERSATION_CONTEXT:
            CONVERSATION_CONTEXT[session_id] = []
        
        CONVERSATION_CONTEXT[session_id].append({"role": "user", "content": mensaje})
        
        mensajes = [{"role": "system", "content": obtener_prompt_sistema(rol)}]
        mensajes.extend(CONVERSATION_CONTEXT[session_id][-10:])
        
        chat_completion = groq_client.chat.completions.create(
            messages=mensajes,
            model="llama-3.1-8b-instant",
            temperature=0.7,
            max_tokens=1500,
            top_p=0.9
        )
        
        respuesta = chat_completion.choices[0].message.content
        CONVERSATION_CONTEXT[session_id].append({"role": "assistant", "content": respuesta})
        
        return jsonify({"response": respuesta, "role": rol, "session_id": session_id})
        
    except Exception as e:
        print(f"Error en chat: {e}")
        return jsonify({"error": str(e)}), 500

@app.route("/api/sugerencias-preguntas")
def api_sugerencias_preguntas():
    """Preguntas sugeridas."""
    return jsonify({"preguntas": [
        "¿Cuál es el significado del contenido GC?",
        "¿Por qué E. coli K-12 es importante?",
        "¿Qué indica la densidad de genes?"
    ]})

@app.route("/")
def index():
    """Página principal con análisis del genoma seleccionado."""
    global CEPA_ACTUAL
    cepa_param = request.args.get('cepa', 'NC_000913.3')
    CEPA_ACTUAL = cepa_param

    registro = descargar_genoma(CEPA_ACTUAL)
    if not registro:
        return "Error: No se pudo descargar el genoma", 500

    secuencia_genoma = str(registro.seq).upper()
    genes = extraer_genes_detallados(registro)
    genes_ordenados = sorted(genes, key=lambda x: x['longitud'], reverse=True)
    compactacion = compactacion_genica(registro)
    espacios = analizar_espacios_basura(registro, genes_ordenados)
    
    # CONTEO DUAL: CDS (correcto) y Genoma Completo (exploratorio)
    conteo_cds = contar_codones_en_cds(registro)
    conteo_genoma = contar_tripletes_genoma_completo(secuencia_genoma)
    stops_finales = contar_stops_finales(registro)
    
    # Contar ATG iniciales (uno por gen)
    atg_iniciales, genes_sin_atg = contar_atg_iniciales(registro)
    
    # Validación con literatura (usando stops_finales)
    validacion = validar_con_literatura(registro, stops_finales)
    
    # Estadísticas de ATG
    atg_stats = {
        "total_en_cds": conteo_cds.get("ATG", 0),  # Todos los ATG dentro de genes
        "iniciales": atg_iniciales,  # Solo los ATG de inicio (uno por gen)
        "no_iniciales": conteo_cds.get("ATG", 0) - atg_iniciales,  # ATG internos
        "total_genoma": conteo_genoma.get("ATG", 0),  # Todos los ATG del genoma
        "genes_sin_atg": len(genes_sin_atg)
    }
    
    # Estadísticas CDS con stops_finales
    codones_stop_cds = stops_finales
    
    total_stop_cds = sum(codones_stop_cds.values())
    proporciones_stop_cds = {
        k: (v / total_stop_cds) * 100 if total_stop_cds > 0 else 0 
        for k, v in codones_stop_cds.items()
    }
    
    # Estadísticas Genoma Completo (exploratorias)
    codones_stop_genoma = {
        "TAA": conteo_genoma.get("TAA", 0),
        "TAG": conteo_genoma.get("TAG", 0),
        "TGA": conteo_genoma.get("TGA", 0)
    }
    
    total_stop_genoma = sum(codones_stop_genoma.values())
    proporciones_stop_genoma = {
        k: (v / total_stop_genoma) * 100 if total_stop_genoma > 0 else 0 
        for k, v in codones_stop_genoma.items()
    }
    
    return render_template(
        "index.html",
        cepa_actual=CEPA_ACTUAL,
        genome_length=len(secuencia_genoma),
        gc=contenido_gc(secuencia_genoma),
        num_genes=len(genes),
        gene_density=compactacion["densidad"],
        compactacion=compactacion,
        gen_mas_grande=genes_ordenados[0] if genes_ordenados else None,
        gen_mas_pequeno=genes_ordenados[-1] if genes_ordenados else None,
        top_genes=genes_ordenados[:15],
        
        # CDS (correcto)
        atg_stats=atg_stats,
        codones_stop_cds=codones_stop_cds,
        stop_props_cds=proporciones_stop_cds,
        stops_finales=stops_finales,
        
        # Genoma completo (exploratorio)
        atg_total_genoma=conteo_genoma.get("ATG", 0),
        codones_stop_genoma=codones_stop_genoma,
        stop_props_genoma=proporciones_stop_genoma,
        
        # Validación
        validacion=validacion,
        
        espacios=espacios[:10],
        genomas_disponibles=list(GENOMES.keys())
    )

@app.route("/api/genes")
def api_genes():
    """API para obtener todos los genes con búsqueda."""
    registro = descargar_genoma(CEPA_ACTUAL)
    if not registro:
        return jsonify({"error": "Error loading genome"}), 500

    genes = extraer_genes_detallados(registro)

    search = request.args.get('search', '').lower()
    if search:
        genes = [g for g in genes if search in g['gen'].lower() or search in g['producto'].lower()]

    # Ordenar por longitud SOLO para display
    genes_ordenados = sorted(genes, key=lambda x: x['longitud'], reverse=True)

    return jsonify({
        "total": len(genes_ordenados),
        "genes": genes_ordenados[:100]
    })

@app.route("/api/gen/<int:gen_index>")
def api_gen_detalle(gen_index):
    """
    API para ver detalles completos de un gen específico.
    Usa index_original para encontrar el gen correcto.
    """
    registro = descargar_genoma(CEPA_ACTUAL)
    if not registro:
        return jsonify({"error": "Error loading genome"}), 500

    genes = extraer_genes_detallados(registro)

    # Buscar gen por index_original
    gen = next((g for g in genes if g['index_original'] == gen_index), None)

    if not gen:
        return jsonify({"error": "Gene not found"}), 404
    secuencia = gen['secuencia']
    
    # Dividir en codones
    codones_gen = [secuencia[i:i+3] for i in range(0, len(secuencia) - 2, 3)]
    contador_codones = Counter(codones_gen)
    
    # Tabla genética estándar
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
    
    # Traducir a proteína
    proteina = ''.join([tabla_genetica.get(codones_gen[i], '?') for i in range(len(codones_gen))])
    
    return jsonify({
        "gen": gen,
        "num_codones": len(codones_gen),
        "codones_frecuentes": dict(contador_codones.most_common(10)),
        "proteina_sintesis": proteina,
        "tamaño_proteina": len(proteina),
        "codones_detalles": list(contador_codones.items())
    })

@app.route("/api/datos-graficos")
def api_datos_graficos():
    """API para obtener datos para gráficos de Chart.js."""
    registro = descargar_genoma(CEPA_ACTUAL)
    if not registro:
        return jsonify({"error": "Error loading genome"}), 500
    
    secuencia_genoma = str(registro.seq).upper()
    conteo_cds = contar_codones_en_cds(registro)
    conteo_genoma = contar_tripletes_genoma_completo(secuencia_genoma)
    stops_finales = contar_stops_finales(registro)
    atg_iniciales, _ = contar_atg_iniciales(registro)
    
    # Datos para gráfico de comparación
    datos_comparacion = {
        "labels": ["ATG Iniciales", "TAA", "TAG", "TGA"],
        "cds": [
            atg_iniciales,
            stops_finales.get("TAA", 0),
            stops_finales.get("TAG", 0),
            stops_finales.get("TGA", 0)
        ],
        "genoma": [
            conteo_genoma.get("ATG", 0),
            conteo_genoma.get("TAA", 0),
            conteo_genoma.get("TAG", 0),
            conteo_genoma.get("TGA", 0)
        ]
    }
    
    # Datos para gráfico de proporciones
    codones_stop_cds = stops_finales
    
    return jsonify({
        "comparacion": datos_comparacion,
        "proporciones_stop": codones_stop_cds
    })

@app.route("/api/cambiar-cepa", methods=['POST'])
def cambiar_cepa():
    """Cambia la cepa actual del análisis."""
    global CEPA_ACTUAL
    data = request.get_json()
    nueva_cepa = data.get('cepa_id')

    if not nueva_cepa:
        return jsonify({"error": "cepa_id requerido"}), 400

    # Validar que la cepa existe
    test = descargar_genoma(nueva_cepa)
    if not test:
        return jsonify({"error": "Cepa no válida o no encontrada"}), 404

    CEPA_ACTUAL = nueva_cepa
    return jsonify({"success": True, "cepa_actual": CEPA_ACTUAL})

@app.route("/api/comparar-genomas")
def api_comparar_genomas():
    """
    API para comparar 2 genomas específicos.
    FUNCIONALIDAD del Proyecto 1.
    """
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
    
    codones1 = contar_codones_en_cds(reg1)
    codones2 = contar_codones_en_cds(reg2)
    
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
    
    diferencias_codones.sort(key=lambda x: abs(x['diferencia']), reverse=True)
    
    comparativa = {
        "genoma1": {
            "nombre": genoma1_nombre,
            "longitud": len(seq1),
            "num_genes": len(genes1),
            "gc": contenido_gc(seq1),
            "compactacion": (sum(g['longitud'] for g in genes1) / len(seq1)) * 100,
        },
        "genoma2": {
            "nombre": genoma2_nombre,
            "longitud": len(seq2),
            "num_genes": len(genes2),
            "gc": contenido_gc(seq2),
            "compactacion": (sum(g['longitud'] for g in genes2) / len(seq2)) * 100,
        },
        "diferencias_codones": diferencias_codones[:30],
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
    """
    Retorna lista de genomas disponibles.
    FUNCIONALIDAD del Proyecto 1.
    """
    return jsonify({
        "genomas": list(GENOMES.keys()),
        "total": len(GENOMES)
    })

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=True)