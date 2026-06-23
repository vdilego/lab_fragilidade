# Roteiro Didático — Laboratório de Heterogeneidade e Fragilidade

**Disciplina**: Métodos Demográficos Avançados  
**Duração**: 2 sessões × 1h30  
**Plataforma**: Posit Cloud  
**Prof. Vanessa Di Lego** · CEDEPLAR/UFMG

---

## Sessão 1 — Heterogeneidade e Modelos de Fragilidade

### Objetivos de aprendizagem

Ao final desta sessão, o/a aluno/a deverá ser capaz de:

1. Simular populações heterogêneas com fragilidade gama e interpretar o papel de σ²
2. Demonstrar numericamente que µ̄(x) desacelera mesmo quando µ₀(x) é Gompertziano crescente
3. Estimar e comparar modelos de fragilidade gama e log-normal via `frailtypack`
4. Quantificar o viés de atenuação no modelo de Cox quando a fragilidade é ignorada
5. Interpretar θ (variância da fragilidade) e relacioná-lo ao plateau de mortalidade (b/σ²)

### Estrutura (90 min)

| Tempo | Atividade |
|-------|-----------|
| 0–15  | Motivação: o problema central (slides 2–5 da aula teórica) · simulação de três populações |
| 15–35 | Visualização da desaceleração de µ̄(x) vs. µ₀(x) · distribuição de Z |
| 35–55 | Cox sem fragilidade → viés · Cox com fragilidade gama (`frailty()`) |
| 55–75 | Comparação gama vs. log-normal: AIC, caudas, identificabilidade |
| 75–90 | Exercício: plateau teórico b/σ² e Z̄(x) por σ² |

### Perguntas de discussão

- Por que a curva observada µ̄(x) pode *declinar* mesmo quando µ₀(x) é monotonamente crescente?
- O plateau de mortalidade em centenários é evidência de desaceleração biológica ou artefato de seleção? Como distinguir?
- Por que o viés de atenuação no Cox é *subestimação* (não superestimação) do HR?
- Em que contexto a fragilidade log-normal seria preferível à gama? O que o AIC diz nos dados simulados?

---

## Sessão 2 — Riscos Competitivos, Fine & Gray e Dados Reais

### Objetivos de aprendizagem

1. Construir a variável de desfecho composto (0=censura, 1=causa k, 2=causa competidora)
2. Calcular e interpretar a CIF de Aalen-Johansen via `survfit(Surv(t, factor(status)))`
3. Demonstrar que 1-KM superestima a incidência cumulativa na presença de riscos competitivos
4. Estimar o modelo de Fine & Gray via `finegray()` e interpretar o SHR
5. Distinguir quando usar causa-específica (etiologia) vs. Fine & Gray (prognóstico)
6. Aplicar fragilidade a dados reais (`pbc`, `veteran`) e interpretar θ e os Z estimados

### Estrutura (90 min)

| Tempo | Atividade |
|-------|-----------|
| 0–15  | `mgus2`: construção do desfecho competitivo · estatísticas descritivas |
| 15–35 | CIF vs. 1-KM: a superestimação · verificação numérica CIF₁+CIF₂+S=1 |
| 35–55 | Fine & Gray via `finegray()` · tabela CS-HR vs. SHR · inversão do efeito da idade |
| 55–70 | `pbc`: fragilidade gama por estágio · θ, Z por estágio, HRs com e sem fragilidade |
| 70–90 | `veteran`: efeito fixo vs. fragilidade · exercício integrador |

### Perguntas de discussão

- O CS-HR da idade para progressão em `mgus2` é 1.016, mas o SHR é 0.983 — sinais opostos. Como isso é possível? Qual o mecanismo?
- No `pbc`, o tratamento com D-penicilamina tem HR ≈ 1.00 em ambos os modelos. O que isso significa historicamente? Como a fragilidade afeta a inferência sobre o tratamento?
- No `veteran`, o tipo celular *smallcell* tem Z = 1.16 e *adeno* Z = 1.45. Por que modelar isso como fragilidade (efeito aleatório) em vez de efeito fixo? Quando a distinção importa?
- Qual a conexão entre o cruzamento de curvas de mortalidade negro/branco nos EUA (slide 5 da aula) e o que observamos nos dados `mgus2` entre sexos?

---

## Exercício Final Integrador

**Contexto** (baseado no slide 44 da aula teórica):

Uma pesquisadora estuda mortalidade por câncer de mama (1970–2010). Após ajustar por estágio e tratamento, a mortalidade de mulheres de baixa renda começa mais alta mas **cruza** a de alta renda por volta dos 75 anos.

1. Quais explicações são possíveis para o cruzamento? Como distingui-las empiricamente?
2. Escreva o modelo de fragilidade univariado e identifique todos os parâmetros.
3. O modelo seria identificável? O que seria necessário para torná-lo identificável?
4. Qual seria o efeito de ignorar a fragilidade nos coeficientes de renda do Cox?
5. Que dado adicional permitiria usar um modelo de fragilidade correlacionada? Que nova pergunta poderíamos responder?

*Dica*: rode a Parte 4 do script `lab1_fragilidade.R`, que simula exatamente este cenário.

---

## Referências

### Leitura obrigatória (antes do lab)

- Vaupel JW, Manton KG, Stallard E (1979). *Demography* **16**(3):439–454.
- Vaupel JW, Yashin AI (1985). *The American Statistician* **39**(3):176–185.
- Fine JP, Gray RJ (1999). *JASA* **94**(446):496–509.

### Leitura complementar

- Vaupel JW, Yashin AI (1987). *Demography* **24**(1):123–135. [Ressuscitação repetida]
- Vaupel JW, Yashin AI, Manton KG (1988). *Mathematical Population Studies* **1**(1):21–48. [Efeito Phoenix]
- Elbers C, Ridder G (1982). *Review of Economic Studies* **49**(3):403–409. [Identificabilidade]
- Wienke A (2003). MPIDR Working Paper WP 2003-032. [Distribuições de fragilidade]
- Duchateau L, Janssen P (2008). *The Frailty Model*. Springer. [Cap. 2, 4 e 7]
- Beyersmann J, Allignol A, Schumacher M (2012). *Competing Risks and Multistate Models with R*. Springer. [Cap. 2 e 4]
