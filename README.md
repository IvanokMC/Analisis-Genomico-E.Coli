# 🧬 Análisis Genómico de *Escherichia coli* K-12 MG1655

Proyecto de Bioinformática para el análisis completo del genoma de *E. coli* K-12 MG1655 utilizando Python, BioPython, Flask y **Asistente IA integrado**.

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue)](https://www.python.org/)
[![Flask](https://img.shields.io/badge/Flask-3.0-green)](https://flask.palletsprojects.com/)
[![BioPython](https://img.shields.io/badge/BioPython-1.83-orange)](https://biopython.org/)
[![AI](https://img.shields.io/badge/AI-Llama%203.1%2070B-purple)](https://groq.com/)
[![License](https://img.shields.io/badge/License-Academic-yellow.svg)](README.md)

## 📋 Descripción del Proyecto

Sistema web interactivo para el análisis genómico de *Escherichia coli* K-12 MG1655 con **Asistente IA integrado** que proporciona interpretación en tiempo real de los datos genómicos.

### 🌟 Características Principales

#### Análisis Genómico Completo
- ✅ **Análisis de codones**: Conteo de codones de inicio (ATG) y terminación (TAA, TAG, TGA)
- ✅ **Validación científica**: Comparación con valores reportados en la literatura (Blattner et al., 1997)
- ✅ **Análisis dual**: Método biológicamente correcto (CDS) vs análisis exploratorio (genoma completo)
- ✅ **Compactación génica**: Cálculo de densidad y porcentaje codificante
- ✅ **Comparativa de genomas**: Análisis comparativo entre diferentes cepas de *E. coli*
- ✅ **Visualización interactiva**: Interfaz web con gráficos y estadísticas detalladas

#### 🤖 **Asistente IA (NUEVO)**
- 🎯 **Chat inteligente** con Llama-3.1-8b-instant vía Groq API
- 🎓 **4 niveles de explicación** adaptativos:
  - **Estudiante**: Explicaciones simples y didácticas con analogías
  - **Investigador**: Análisis técnico profundo con terminología científica
  - **Docente**: Balance entre rigor científico y claridad pedagógica
  - **Divulgador**: Lenguaje accesible para audiencia general
- 💡 **Sugerencias de preguntas** categorizadas
- 📊 **Contexto genómico en tiempo real** - La IA conoce los datos actuales del análisis
- 🧠 **Memoria de conversación** - Mantiene contexto entre preguntas
- 🔬 **Interpretación de resultados** - Análisis contextual de comparaciones y genes

## 🎯 Objetivos

1. Identificar y cuantificar codones de inicio y terminación en el genoma
2. Calcular la densidad génica y porcentaje de compactación
3. Validar resultados contra valores de referencia en la literatura científica
4. Comparar diferentes métodos de conteo de codones
5. Analizar espacios intergénicos (regiones no codificantes)
6. Comparar múltiples cepas de *E. coli*
7. **Proporcionar interpretación asistida por IA** de los resultados genómicos

## 🔬 Metodología

### Método Biológicamente Correcto (CDS)

Conteo de codones **únicamente dentro de secuencias codificantes (CDS) anotadas**, respetando el marco de lectura correcto:

- ✅ Solo analiza genes anotados por expertos
- ✅ Respeta el marco de lectura (+0)
- ✅ Los resultados coinciden con la literatura científica

### Método Exploratorio (Genoma Completo)

Conteo de **todos los tripletes** en el genoma completo:

- 📊 Incluye regiones intergénicas
- 📊 Incluye todos los marcos de lectura (+0, +1, +2)
- 📊 Útil para análisis estadístico comparativo

## 📊 Resultados Esperados (K-12 MG1655)

- **Longitud**: ~4.64 Mb
- **Contenido GC**: ~50.79%
- **Genes**: ~4,318
- **Compactación**: ~86.76%
- **Codones STOP en CDS**:
  - TGA: ~61%
  - TAA: ~30%
  - TAG: ~9%

## 🚀 Instalación

### Requisitos

- Python 3.10 o superior
- Sistema operativo: Linux (Ubuntu 22.04+), macOS o Windows (WSL2)
- Conexión a internet

### Pasos

```bash
# 1. Clonar el repositorio
git clone https://github.com/[tu-usuario]/analisis-genoma-ecoli.git
cd analisis-genoma-ecoli

# 2. Crear y activar entorno virtual
python3 -m venv venv
source venv/bin/activate    # Linux/macOS
# venv\Scripts\activate     # Windows

# 3. Instalar dependencias
pip install -r requirements.txt

# 4. Ejecutar la aplicación
python app.py
```

Accede en: **http://localhost:5000**

## 🛠️ Tecnologías

- **Backend**: Flask, BioPython, Groq API
- **Frontend**: HTML, CSS, JavaScript, Chart.js
- **IA**: Llama-3.1-8b-instant (vía Groq)
- **Datos**: NCBI GenBank

## 📡 API Endpoints

### Análisis Genómico
| Ruta | Método | Descripción |
|------|--------|-------------|
| `/` | GET | Página principal con dashboard |
| `/api/genes` | GET | Lista de genes (con `?search=`) |
| `/api/gen/<id>` | GET | Detalles de un gen específico |
| `/api/datos-graficos` | GET | Datos para gráficos |
| `/api/comparar-genomas` | GET | Comparativa entre cepas |
| `/api/genomas-disponibles` | GET | Lista de genomas |

### 🤖 Asistente IA
| Ruta | Método | Descripción |
|------|--------|-------------|
| `/api/chat` | POST | Enviar mensaje al chatbot |
| `/api/sugerencias-preguntas` | GET | Preguntas sugeridas |

## 🤖 Uso del Asistente IA

### 1. Seleccionar Nivel
- **🎓 Estudiante**: Explicaciones simples
- **🔬 Investigador**: Análisis técnico
- **👨‍🏫 Docente**: Enfoque pedagógico
- **📢 Divulgador**: Lenguaje accesible

### 2. Hacer Preguntas

**Ejemplos:**
- "¿Qué es el contenido GC?"
- "Resume los datos del genoma"
- "¿Cómo optimizar un gen para expresión?"
- "¿Qué diferencias hay entre K-12 y O157:H7?"

### 3. Usar Sugerencias
Haz clic en las preguntas sugeridas para enviarlas automáticamente.

## 🧪 Validación

Resultados validados contra **Blattner et al. (1997)**:

| Parámetro | Esperado | Obtenido | Estado |
|-----------|----------|----------|--------|
| Longitud | 4.5-4.7 Mb | ~4.64 Mb | ✅ |
| GC | 50.5-51.0% | ~50.79% | ✅ |
| Genes | 4,200-4,400 | ~4,318 | ✅ |
| TAA | ~30% | ~29.78% | ✅ |
| TAG | ~9% | ~9.22% | ✅ |
| TGA | ~61% | ~61.00% | ✅ |

## 📚 Referencias

1. Blattner, F. R., et al. (1997). "The complete genome sequence of Escherichia coli K-12." *Science*, 277(5331), 1453-1462.
2. Riley, M., et al. (2006). "Escherichia coli K-12: a cooperatively developed annotation snapshot—2005." *Nucleic Acids Research*, 34(1), 1-9.
3. Keseler, I. M., et al. (2017). "The EcoCyc database." *Nucleic Acids Research*, 45(D1), D543-D550.

## 🎓 Interpretación

### Diferencia entre Métodos
- Solo **5.9%** de ATG en genoma son inicios reales
- Solo **2.4%** de STOP en genoma son finales reales
- Demuestra la importancia de usar anotaciones CDS

### Compactación
- ~86.76% indica alta eficiencia evolutiva
- Típico de genomas procariotas
- Poco "ADN basura"

## ⚠️ Notas

- **Primera carga**: ~10-30 segundos (descarga de NCBI)
- **Caché**: Genomas se guardan en memoria
- **Email NCBI**: Cambiar en `app.py`
- **API Key**: Ya incluida para desarrollo

## 🛠️ Solución de Problemas

### "No module named 'Bio'"
```bash
source venv/bin/activate
pip install biopython
```

### "Address already in use"
```bash
kill $(lsof -t -i:5000)
```

### El chat no responde
- Verificar conexión a internet
- Revisar consola del navegador (F12)
- Verificar logs del servidor

## 👨‍💻 Autores

- Diego Shaid Ninancuro Huarayo
- Luis Ángel Mogrovejo Ccorimanya
- George Ivanok Muñoz Castillo
- Gustavo Pantoja Olave

**Universidad Nacional San Antonio Abad del Cusco**  
Ingeniería de Sistemas y Computación  
Curso: Bioinformática

## 📄 Licencia

Proyecto académico – Uso educativo y no comercial.

## 🙏 Agradecimientos

- **NCBI** por acceso público a genomas
- **BioPython** por herramientas de análisis
- **Flask** por el framework web
- **Chart.js** por visualización
- **Groq** por acceso a llama-3.1-8b-instant

---

<div align="center">

**🧬 Desarrollado con ❤️ para el avance de la bioinformática**

[![UNSAAC](https://img.shields.io/badge/UNSAAC-Ingeniería%20de%20Sistemas-blue)](https://www.unsaac.edu.pe/)

**Última actualización: Febrero 2026**

</div>