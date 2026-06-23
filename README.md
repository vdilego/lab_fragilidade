# Laboratório: Heterogeneidade, Fragilidade e Riscos Competitivos

**Disciplina**: Métodos Demográficos Avançados — Doutorado  
**Repositório**: [github.com/vdilego/lab_fragilidade](https://github.com/vdilego/lab_fragilidade)  
**Plataforma**: [Posit Cloud](https://posit.cloud)  
**Prof. Vanessa Di Lego** · CEDEPLAR/UFMG · Junho 2026

---

## Estrutura do Laboratório

O laboratório está dividido em **duas sessões de 1h30** e acompanha a aula teórica
*Heterogeneidade das Idades à Morte e Modelos de Fragilidade* (CEDEPLAR, 16/06/2026).

| Sessão | Tema | Script |
|--------|------|--------|
| Lab 1 | Heterogeneidade e modelos de fragilidade com dados simulados | `lab1_fragilidade.R` |
| Lab 2 | Riscos competitivos, Fine & Gray e dados reais | `lab2_dados_reais.R` |

Scripts auxiliares:

- `00_setup.R` — instalação de pacotes e tema global ggplot

---

## Pacotes necessários

```r
install.packages(c(
  "survival",      # Cox, fragilidade, CIF, Fine & Gray (finegray)
  "frailtypack",   # fragilidade gama/log-normal semiparamétrica
  "flexsurv",      # modelos paramétricos com fragilidade
  "cmprsk",        # Fine & Gray clássico (crr)
  "tidyverse",
  "ggplot2",
  "patchwork",
  "broom",
  "ggsci",
  "scales"
))
```

> **Nota**: `finegray()` e `survfit(Surv(t, factor(status)))` para CIF Aalen-Johansen
> estão disponíveis no pacote `survival` a partir da versão 3.0 (CRAN, sem dependências extras).

---

## Dados

Todos os datasets usados são **embutidos no pacote `survival`** — nenhum arquivo externo é necessário.

| Dataset | Fonte | N | Uso no laboratório |
|---------|-------|---|--------------------|
| `mgus2` | Kyle et al. (2002) *NEJM* | 1.384 | CIF, Fine & Gray — progressão vs. morte |
| `pbc`   | Fleming & Harrington (1991) | 312 | Fragilidade gama por estágio clínico |
| `veteran` | Kalbfleisch & Prentice (1980) | 137 | Fragilidade por tipo celular |

---

## Slides

Os slides com síntese dos resultados estão disponíveis em:

- [`slides_resultados_lab.html`](slides_resultados_lab.html) — apresentação Reveal.js (abrir no navegador)

Para visualizar via GitHub Pages, acesse:  
`https://vdilego.github.io/lab_fragilidade/slides_resultados_lab.html`

---

## Roteiro didático

O roteiro completo com objetivos de aprendizagem, estrutura das sessões e
perguntas de discussão está em [`roteiro_didatico.md`](roteiro_didatico.md).

---

## Referências bibliográficas

### Modelos de fragilidade

- Vaupel JW, Manton KG, Stallard E (1979). The impact of heterogeneity in individual frailty on the dynamics of mortality. *Demography* **16**(3):439–454.
- Vaupel JW, Yashin AI (1985). Heterogeneity's ruses: some surprising effects of selection on population dynamics. *The American Statistician* **39**(3):176–185.
- Vaupel JW, Yashin AI (1987). Repeated resuscitation: how lifesaving alters life tables. *Demography* **24**(1):123–135.
- Vaupel JW, Yashin AI (1987). Targeting lifesaving: demographic linkages between population structure and life expectancy. *European Journal of Population* **2**(3):335–360.
- Vaupel JW, Yashin AI, Manton KG (1988). Debilitation's aftermath: stochastic process models of mortality. *Mathematical Population Studies* **1**(1):21–48.
- Elbers C, Ridder G (1982). True and spurious duration dependence. *Review of Economic Studies* **49**(3):403–409.
- Wienke A (2003). Frailty models. *MPIDR Working Paper* WP 2003-032.
- Duchateau L, Janssen P (2008). *The Frailty Model*. Springer.

### Riscos competitivos e Fine & Gray

- Fine JP, Gray RJ (1999). A proportional hazards model for the subdistribution of a competing risk. *JASA* **94**(446):496–509. DOI: [10.1080/01621459.1999.10474144](https://doi.org/10.1080/01621459.1999.10474144)
- Beyersmann J, Allignol A, Schumacher M (2012). *Competing Risks and Multistate Models with R*. Springer.

### Dados

- Kyle RA et al. (2002). A long-term study of prognosis in monoclonal gammopathy of undetermined significance. *NEJM* **346**(8):564–569.
- Fleming TR, Harrington DP (1991). *Counting Processes and Survival Analysis*. Wiley.
- Kalbfleisch JD, Prentice RL (1980). *The Statistical Analysis of Failure Time Data*. Wiley.
