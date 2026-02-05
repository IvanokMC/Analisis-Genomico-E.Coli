# 🧬 Análisis Genómico de *Escherichia coli* K-12 MG1655

Proyecto de Bioinformática para análisis completo del genoma de *E. coli* K-12 MG1655 utilizando Python, BioPython y Flask.

## 📋 Descripción del Proyecto

Este proyecto implementa un análisis genómico completo de *Escherichia coli* K-12 MG1655, incluyendo:

- **Análisis de codones**: Conteo de codones de inicio (ATG) y terminación (TAA, TAG, TGA)
- **Validación científica**: Comparación con valores reportados en la literatura (Blattner et al., 1997)
- **Análisis dual**: Método biológicamente correcto (CDS) vs análisis exploratorio (genoma completo)
- **Compactación génica**: Cálculo de densidad y porcentaje codificante
- **Comparativa de genomas**: Análisis comparativo entre diferentes cepas de *E. coli*
- **Visualización interactiva**: Interfaz web con gráficos y estadísticas detalladas

## 🎯 Objetivos

1. Identificar y cuantificar codones de inicio y terminación en el genoma
2. Calcular la densidad génica y porcentaje de compactación
3. Validar resultados contra valores de referencia en la literatura científica
4. Comparar diferentes métodos de conteo de codones
5. Analizar espacios intergénicos (regiones no codificantes)
6. Comparar múltiples cepas de *E. coli*

## 🔬 Metodología

### Método Biológicamente Correcto (CDS)

Conteo de codones **únicamente dentro de secuencias codificantes (CDS) anotadas**, respetando el marco de lectura correcto:

- Solo analiza genes anotados por expertos
- Respeta el marco de lectura (+0)
- Los resultados coinciden con la literatura científica

### Método Exploratorio (Genoma Completo)

Conteo de **todos los tripletes** en el genoma completo:

- Incluye regiones intergénicas
- Incluye todos los marcos de lectura (+0, +1, +2)
- Útil para análisis estadístico comparativo

## 📊 Resultados Esperados

### Información General
- **Longitud del genoma**: ~4.64 Mb
- **Contenido GC**: ~50.79%
- **Total de genes**: ~4,318
- **Densidad génica**: ~930 genes/Mb
- **Compactación**: ~86.76%

### Codones de Terminación (CDS)
- **TAA**: ~1,286 (29.78%)
- **TAG**: ~398 (9.22%)
- **TGA**: ~2,634 (61.00%)

Estos valores coinciden con los reportados por Blattner et al. (1997).

## 🚀 Instalación

### Requisitos Previos

- Ubuntu 22.04 LTS (o superior)
- Python 3.10+
- Conexión a internet (para descargar genomas de NCBI)

### Instalación Paso a Paso
```bash
# 1. Actualizar sistema
sudo apt update && sudo apt upgrade -y

# 2. Instalar dependencias del sistema
sudo apt install python3 python3-pip python3-venv python3-dev build-essential git -y

# 3. Clonar repositorio (si aplica)
git clone https://github.com/tu-usuario/proyecto-genoma-ecoli.git
cd proyecto-genoma-ecoli

# 4. Crear entorno virtual
python3 -m venv venv

# 5. Activar entorno virtual
source venv/bin/activate

# 6. Instalar dependencias de Python
pip install -r requirements.txt
```

## ▶️ Ejecución

### Modo Desarrollo
```bash
# Activar entorno virtual
source venv/bin/activate

# Ejecutar aplicación
python app.py
```

### Modo Producción (AWS EC2)
```bash
# Activar entorno virtual
source venv/bin/activate

# Ejecutar con nohup para mantener en background
nohup python app.py > app.log 2>&1 &

# Ver el log
tail -f app.log
```

### Detener la Aplicación
```bash
# Encontrar el proceso
ps aux | grep app.py

# Detener (reemplaza PID con el número del proceso)
kill PID
```

## 📦 Dependencias
```
Flask==3.0.0              # Framework web
biopython==1.83           # Análisis bioinformático
pandas==2.1.4             # Manejo de datos
matplotlib==3.8.2         # Gráficos (opcional)
```

## 🔧 Funciones Principales

### `descargar_genoma(genome_id)`
Descarga un genoma desde NCBI GenBank usando BioPython.

### `contar_codones_en_cds(registro)`
Cuenta codones solo dentro de CDS anotados (método correcto).

### `contar_tripletes_genoma_completo(secuencia)`
Cuenta todos los tripletes en el genoma (método exploratorio).

### `validar_con_literatura(registro, conteo_cds)`
Valida resultados contra valores esperados de Blattner et al. (1997).

### `compactacion_genica(registro)`
Calcula porcentaje de compactación y densidad génica.

### `extraer_genes_detallados(registro)`
Extrae información completa de todos los genes anotados.

## 📈 API Endpoints

| Endpoint | Método | Descripción |
|----------|--------|-------------|
| `/` | GET | Página principal con análisis completo |
| `/api/genes` | GET | Lista de genes (con búsqueda opcional) |
| `/api/datos-graficos` | GET | Datos para gráficos de Chart.js |
| `/api/comparar-genomas` | GET | Comparación entre 2 genomas |
| `/api/genomas-disponibles` | GET | Lista de genomas disponibles |

### Ejemplos de Uso
```bash
# Obtener todos los genes
curl http://localhost:5000/api/genes

# Buscar genes específicos
curl http://localhost:5000/api/genes?search=lac

# Comparar genomas
curl "http://localhost:5000/api/comparar-genomas?genoma1=K-12%20MG1655%20(Referencia)&genoma2=K-12%20W3110"

# Obtener datos para gráficos
curl http://localhost:5000/api/datos-graficos
```

## 🧪 Validación

Los resultados fueron validados contra:

**Referencia Principal:**
> Blattner, F. R., et al. (1997). "The complete genome sequence of Escherichia coli K-12." *Science*, 277(5331), 1453-1462.

**Valores Esperados:**
- Longitud: 4.5-4.7 Mb ✓
- Contenido GC: 50.5-51.0% ✓
- Número de genes: 4,200-4,400 ✓
- Proporción TAA: ~30% ✓
- Proporción TAG: ~9% ✓
- Proporción TGA: ~61% ✓

## 📚 Referencias

1. Blattner, F. R., et al. (1997). "The complete genome sequence of Escherichia coli K-12." *Science*, 277(5331), 1453-1462.

2. Riley, M., et al. (2006). "Escherichia coli K-12: a cooperatively developed annotation snapshot—2005." *Nucleic Acids Research*, 34(1), 1-9.

3. Keseler, I. M., et al. (2017). "The EcoCyc database: reflecting new knowledge about Escherichia coli K-12." *Nucleic Acids Research*, 45(D1), D543-D550.

## 🎓 Interpretación de Resultados

### Diferencia entre Métodos de Conteo

El análisis revela que:

- **Solo 5.9%** de los tripletes ATG en el genoma son inicios reales de genes
- **Solo 2.4%** de los tripletes STOP son finales reales de genes
- El resto aparece por **distribución estadística** en regiones no codificantes

Esta diferencia demuestra la importancia de usar anotaciones expertas (CDS) para análisis biológicamente relevantes.

### Compactación Génica

*E. coli* presenta ~86.76% de compactación, lo cual es característico de genomas bacterianos:

- Genomas compactos = alta eficiencia evolutiva
- Poco "ADN basura" = mayoría del genoma codifica proteínas
- Típico de organismos procariotas

### Uso de Codones STOP

La preferencia por TGA (~61%) sobre TAA (~30%) y TAG (~9%) refleja:

- Sesgo evolutivo en el uso de codones
- Optimización para la maquinaria de traducción de *E. coli*
- Consistente con otros estudios de uso de codones

## ⚠️ Notas Importantes

1. **Tiempo de Carga**: La primera carga puede tardar ~10 segundos mientras descarga el genoma de NCBI

2. **Caché**: Los genomas se cachean en memoria para cargas posteriores más rápidas

3. **Seguridad AWS**: Asegúrate de abrir el puerto 5000 en el Security Group de tu instancia EC2

4. **Email NCBI**: El código usa un email de ejemplo. Cámbialo por tu email real en `app.py`:
```python
   Entrez.email = "tu-email@ejemplo.com"
```

## 🐛 Solución de Problemas

### Error: "No module named 'Bio'"
```bash
# Asegúrate de tener el venv activado
source venv/bin/activate
pip install biopython
```

### Error: "Address already in use"
```bash
# El puerto 5000 ya está en uso
# Opción 1: Detener el proceso existente
kill $(lsof -t -i:5000)

# Opción 2: Cambiar el puerto en app.py
app.run(host='0.0.0.0', port=8000, debug=True)
```

### Error: "Connection refused" desde navegador
```bash
# Verifica que el Security Group en AWS tenga el puerto abierto
# EC2 → Security Groups → Inbound Rules → Add Rule
# Type: Custom TCP, Port: 5000, Source: 0.0.0.0/0
```

### Error al descargar genomas
```bash
# Verifica conexión a internet
ping 8.8.8.8

# Verifica que NCBI esté accesible
curl https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi
```


## 👨‍💻 Autores

Luis Ángel Mogrovejo Ccorimanya
Universidad Nacional San Antonio Abad del Cusco
Ingeniería de Sistemas y Computación

George Ivanok Muñoz Castillo
Universidad Nacional San Antonio Abad del Cusco
Ingeniería de Sistemas y Computación

Diego Shaid Ninancuro Huarayo
Universidad Nacional San Antonio Abad del Cusco
Ingeniería de Sistemas y Computación

Wendel Yovan Niño de Guzman Conde
Universidad Nacional San Antonio Abad del Cusco
Ingeniería de Sistemas y Computación

Gustavo Pantoja Olave
Universidad Nacional San Antonio Abad del Cusco
Ingeniería de Sistemas y Computación

## 📄 Licencia

Este proyecto es parte de un trabajo académico para el curso de Bioinformática.

## 🙏 Agradecimientos

- **NCBI** por proporcionar acceso público a genomas
- **BioPython** por las herramientas de análisis bioinformático
- **Flask** por el framework web
- **Chart.js** por las bibliotecas de visualización

---

**Última actualización**: Febrero 2026
