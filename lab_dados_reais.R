################################################################################
##                                                                            ##
##   LABORATÓRIO — PARTE 5: DADOS REAIS                                      ##
##   Heterogeneidade, Fragilidade e Riscos Competitivos                       ##
##                                                                            ##
##   Todos os datasets usados aqui estão embutidos no pacote {survival}      ##
##   — nenhum arquivo externo é necessário.                                   ##
##                                                                            ##
##   Datasets:                                                                ##
##     mgus2    — Gamopatia Monoclonal (Mayo Clinic, n=1384)                  ##
##                Riscos competitivos: progressão vs. morte                   ##
##     pbc      — Cirrose Biliar Primária (Mayo Clinic, n=418)                ##
##                Fragilidade por estágio clínico                             ##
##     veteran  — Câncer de Pulmão — VA Trial (n=137)                        ##
##                Fragilidade por tipo celular                                ##
##                                                                            ##
##   Prof. Vanessa Di Lego — CEDEPLAR/UFMG                                   ##
##                                                                            ##
################################################################################

library(survival)
library(ggplot2)
library(patchwork)
library(tidyverse)

theme_demog <- theme_bw(base_size = 13) +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill = "#2C3E50", color = NA),
        strip.text       = element_text(color = "white", face = "bold"))
theme_set(theme_demog)


## ============================================================================
## DATASET 1 — mgus2: GAMOPATIA MONOCLONAL (MGUS)
## ============================================================================
##
## Contexto clínico e demográfico:
##   MGUS (Monoclonal Gammopathy of Undetermined Significance) é uma condição
##   pré-maligna assintomática. Pacientes correm dois riscos distintos:
##
##   Risco 1 — Progressão para malignidade:  PCM (mieloma múltiplo),
##             amiloidose ou outra doença linfoide relacionada
##   Risco 2 — Morte sem ter progredido (por qualquer outra causa)
##
##   Esses dois desfechos COMPETEM entre si:
##   um paciente que morre não pode mais progredir para PCM (e vice-versa).
##   Ignorar essa competição leva a superestimação da incidência cumulativa.
##
##   Fonte: Kyle RA et al. (2002) NEJM 346(8):564-569.
##          Disponível em survival::mgus2 (n=1384, Mayo Clinic, 1960-1994).

mgus2   <- survival::mgus2 


cat("╔══════════════════════════════════════════════════════╗\n")
cat("║  DATASET: mgus2 — Gamopatia Monoclonal (Mayo Clinic) ║\n")
cat("╚══════════════════════════════════════════════════════╝\n\n")
cat("Variáveis principais:\n")
cat("  id, age, sex, dxyr    → identificação e demográficas\n")
cat("  hgb, creat, mspike    → biomarcadores clínicos\n")
cat("  ptime, pstat          → tempo e status de progressão\n")
cat("  futime, death         → tempo de seguimento e morte\n\n")

cat("Distribuição por sexo:\n")
print(table(mgus2$sex))
cat("\nIdade ao diagnóstico:\n")
print(summary(mgus2$age))
cat("\nM-spike (g/dL):\n")
print(summary(mgus2$mspike))


## ----------------------------------------------------------------------------
## 5.1  Construindo o desfecho composto de risco competitivo
## ----------------------------------------------------------------------------
##
## Convenção usada no vignette oficial do pacote survival:
##   etime = tempo até o PRIMEIRO evento (progressão OU morte), em meses
##   event: 0 = censura administrativa
##          1 = progressão para PCM/amiloidose
##          2 = morte sem progressão (risco competidor)

mgus2 <- mgus2 %>%
  mutate(
    etime = ifelse(pstat == 0, futime, ptime),
    event = ifelse(pstat == 0, 2 * death, 1)
    # Se não progrediu (pstat==0):
    #   event = 2*death → 2 se morreu, 0 se saiu vivo (censura)
    # Se progrediu (pstat==1):
    #   event = 1 (progressão, independente do que aconteceu depois)
  )

cat("\n=== Distribuição do desfecho composto ===\n")
tab_ev <- table(mgus2$event)
names(tab_ev) <- c("0 = censura", "1 = progressão", "2 = morte sem progressão")
print(tab_ev)

cat("\nProporção de progressão em 10 anos (estimativa bruta):\n")
cat(" M:", round(mean(mgus2$event[mgus2$sex == "M"] == 1), 3), "\n")
cat(" F:", round(mean(mgus2$event[mgus2$sex == "F"] == 1), 3), "\n")


## ----------------------------------------------------------------------------
## 5.2  CIF correta vs. 1-KM: a armadilha mais comum
## ----------------------------------------------------------------------------
##
## TEORIA (slide da aula):
##   F_k(t) = ∫₀ᵗ S(u⁻) · h_k(u) du   [CIF de Aalen-Johansen]
##
##   1 - KM_k(t) NÃO é F_k(t) quando há riscos competitivos.
##   O KM trata as mortes por outras causas como CENSURA ALEATÓRIA,
##   o que viola a suposição de independência e SUPERESTIMA F_k(t).
##
## No survival >= 3.0, survfit() com factor(event) calcula CIFs corretas.

## CIF via Aalen-Johansen (survival nativo — CORRETO)
cif_mgus <- survfit(Surv(etime, factor(event)) ~ sex, data = mgus2)

## 1-KM por causa (ERRADO como estimador de incidência)
km_prog_sex <- survfit(Surv(etime, event == 1) ~ sex, data = mgus2)

## Extrair para plot
extrair_cif_survival <- function(obj, n_strata = 2, nomes_strata = c("F", "M")) {
  # Extrai pstate (probabilidade de cada estado) do survfit multi-estado
  # pstate[,1]=sobrevivência, pstate[,2]=causa1, pstate[,3]=causa2

  strata_idx <- obj$strata
  limites    <- c(0, cumsum(strata_idx))

  map_dfr(seq_along(nomes_strata), function(i) {
    idx <- (limites[i] + 1):limites[i + 1]
    tibble(
      time         = obj$time[idx] / 12,   # meses → anos
      cif_prog     = obj$pstate[idx, 2],
      cif_morte    = obj$pstate[idx, 3],
      sobrevivencia = obj$pstate[idx, 1],
      sexo         = nomes_strata[i]
    )
  })
}

df_cif <- extrair_cif_survival(cif_mgus)

## Extrair 1-KM para comparação
df_km <- broom::tidy(km_prog_sex) %>%
  mutate(
    incid_km = 1 - estimate,
    sexo     = ifelse(grepl("F", strata), "F", "M"),
    time     = time / 12
  )

## --- Figura 5.1: 1-KM vs. CIF para progressão ---
fig5_1 <- ggplot() +
  geom_step(
    data = df_km,
    aes(x = time, y = incid_km, color = sexo, linetype = "1-KM (ERRADO)"),
    linewidth = 0.85
  ) +
  geom_step(
    data = df_cif,
    aes(x = time, y = cif_prog, color = sexo, linetype = "CIF Aalen-Johansen (CORRETO)"),
    linewidth = 0.95
  ) +
  scale_color_manual(values = c(F = "#E64B35", M = "#3C5488"),
                     labels = c(F = "Feminino", M = "Masculino"),
                     name   = "Sexo") +
  scale_linetype_manual(
    values = c("1-KM (ERRADO)" = "dashed",
               "CIF Aalen-Johansen (CORRETO)" = "solid"),
    name = "Estimador"
  ) +
  scale_x_continuous(breaks = seq(0, 35, 5)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Figura 5.1 — Progressão para PCM: 1-KM superestima a incidência",
    subtitle = "Morte sem progressão é risco competitivo — não pode ser tratada como censura",
    x = "Anos desde diagnóstico de MGUS",
    y = "Incidência cumulativa de progressão",
    caption  = paste(
      "Dados: mgus2 (n=1384, Mayo Clinic).",
      "Kyle et al. (2002) NEJM 346(8):564-569.",
      "\n1-KM trata mortes como censura aleatória → viola independência → superestima F_k(t)."
    )
  )

fig5_1

## --- Figura 5.2: CIF completa — dois desfechos empilhados ---
fig5_2 <- df_cif %>%
  select(time, sexo, cif_prog, cif_morte) %>%
  pivot_longer(
    cols      = c(cif_prog, cif_morte),
    names_to  = "desfecho",
    values_to = "cif"
  ) %>%
  mutate(
    desfecho = recode(desfecho,
                      cif_prog  = "Progressão para PCM/amiloidose",
                      cif_morte = "Morte sem progressão")
  ) %>%
  ggplot(aes(x = time, y = cif, color = sexo, linetype = desfecho)) +
  geom_step(linewidth = 0.9) +
  scale_color_manual(values = c(F = "#E64B35", M = "#3C5488"),
                     labels = c(F = "Feminino", M = "Masculino"),
                     name   = "Sexo") +
  scale_linetype_manual(
    values = c("Progressão para PCM/amiloidose" = "solid",
               "Morte sem progressão"           = "dashed"),
    name   = "Desfecho"
  ) +
  scale_x_continuous(breaks = seq(0, 35, 5)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Figura 5.2 — CIF por causa e sexo (Aalen-Johansen)",
    subtitle = "Homens têm maior mortalidade geral; mulheres têm maior progressão proporcional",
    x        = "Anos desde diagnóstico",
    y        = "Incidência cumulativa",
    caption  = "CIF₁(t) + CIF₂(t) + S(t) = 1 para cada grupo. Estimativa via survfit(factor(event))."
  )

fig5_2
print(fig5_1 / fig5_2)

## Verificação numérica: soma deve ser 1
cat("\n=== Verificação: CIF₁ + CIF₂ + S(t) = 1? ===\n")
df_cif %>%
  filter(abs(time - round(time)) < 0.05) %>%   # valores anuais
  mutate(soma = cif_prog + cif_morte + sobrevivencia) %>%
  filter(time %in% c(1, 5, 10, 15, 20)) %>%
  select(time, sexo, cif_prog, cif_morte, sobrevivencia, soma) %>%
  mutate(across(where(is.numeric), ~round(.x, 4))) %>%
  print()

## CIF em momentos-chave
cat("\n=== CIF em 5 e 10 anos por sexo ===\n")
s_5  <- summary(cif_mgus, times = 60)
s_10 <- summary(cif_mgus, times = 120)

cat("\n5 anos — Progressão:\n")
print(setNames(round(s_5$pstate[, 2], 3), c("Feminino", "Masculino")))
cat("5 anos — Morte sem progressão:\n")
print(setNames(round(s_5$pstate[, 3], 3), c("Feminino", "Masculino")))

cat("\n10 anos — Progressão:\n")
print(setNames(round(s_10$pstate[, 2], 3), c("Feminino", "Masculino")))
cat("10 anos — Morte sem progressão:\n")
print(setNames(round(s_10$pstate[, 3], 3), c("Feminino", "Masculino")))


## ----------------------------------------------------------------------------
## 5.3  Modelos: Causa-Específica vs. Fine & Gray
## ----------------------------------------------------------------------------
##
## Pergunta substantiva:
##   Quais fatores clínicos aumentam o risco de progressão para PCM?
##   O M-spike elevado — biomarcador da clona maligna — deve ser o principal.
##   Mas como o efeito muda entre CS-Cox e Fine & Gray?

## --- Modelos de causa-específica (Cox por causa) ---
cs_prog <- coxph(Surv(etime, event == 1) ~ age + sex + mspike,
                  data = mgus2)
cs_mort <- coxph(Surv(etime, event == 2) ~ age + sex + mspike,
                  data = mgus2)

## --- Fine & Gray via finegray() nativo do survival ---
##
## finegray() constrói o "weighted pseudo-dataset" onde:
##   - Indivíduos que tiveram o risco competidor permanecem "em risco"
##     para o evento de interesse, mas com peso decrescente
##   - Isso captura o efeito sobre a CIF diretamente
##
## Referência: Fine JP, Gray RJ (1999). A proportional hazards model
##   for the subdistribution of a competing risk. JASA 94(446):496-509.
##   DOI: 10.1080/01621459.1999.10474144

fg_data1 <- finegray(Surv(etime, factor(event)) ~ ., data = mgus2, etype = 1)
fg_prog  <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex + mspike,
                   weight = fgwt, data = fg_data1)

fg_data2 <- finegray(Surv(etime, factor(event)) ~ ., data = mgus2, etype = 2)
fg_mort  <- coxph(Surv(fgstart, fgstop, fgstatus) ~ age + sex + mspike,
                   weight = fgwt, data = fg_data2)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║    TABELA 5.1 — Causa-Específica vs. Fine & Gray (mgus2)            ║\n")
cat("╠══════════════════════════════════════════════════════════════════════╣\n")
cat("║  HR = Hazard Ratio   SHR = Subdistribution Hazard Ratio             ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

tab_comp <- tibble(
  Covariavel          = c("Idade (por ano)", "Sexo masculino", "M-spike (por g/dL)"),
  `HR-CS Progressão`  = round(exp(coef(cs_prog)), 3),
  `SHR Progressão`    = round(exp(coef(fg_prog)), 3),
  `HR-CS Morte`       = round(exp(coef(cs_mort)), 3),
  `SHR Morte`         = round(exp(coef(fg_mort)), 3)
)
print(tab_comp, n = Inf)

cat("\n--- Interpretação linha por linha ---\n\n")

cat("IDADE:\n")
cat("  CS-HR progressão = 1.016  → dado que vivo, cada ano extra de idade\n")
cat("                             aumenta 1.6% o risco instantâneo de progredir\n")
cat("  SHR progressão   = 0.983  → cada ano extra de idade REDUZ 1.7% a CIF\n")
cat("                             Por quê? Mais velhos morrem mais (risco 2),\n")
cat("                             então sobram menos pessoas para progredir.\n")
cat("                             → Efeito indireto via competição!\n\n")

cat("M-SPIKE:\n")
cat("  CS-HR progressão = 2.421  → M-spike dobra aprox. o risco instantâneo\n")
cat("  SHR progressão   = 2.431  → efeito quase idêntico sobre CIF\n")
cat("                             Por quê? M-spike não afeta muito a mortalidade\n")
cat("                             geral, logo a competição não distorce.\n\n")

cat("LIÇÃO CENTRAL:\n")
cat("  Quando a covariável afeta AMBAS as causas, CS-HR ≠ SHR.\n")
cat("  Quando a covariável afeta apenas uma causa, CS-HR ≈ SHR.\n")
cat("  O SHR é o estimador correto para PROGNÓSTICO e COMUNICAÇÃO de risco.\n")
cat("  O CS-HR é o estimador correto para ETIOLOGIA e MECANISMO biológico.\n")


## --- Figura 5.3: CIF ajustada por Fine & Gray ---
## Comparar CIF prevista para perfil médio vs. alto M-spike

perfil_base  <- data.frame(age = 70, sex = factor("F", levels=c("F","M")), mspike = 1.0)
perfil_alto  <- data.frame(age = 70, sex = factor("F", levels=c("F","M")), mspike = 3.0)

# CIF ajustada: predita a partir do Fine & Gray (via predict no pseudo-dataset)
# Alternativa mais direta: refit do survfit no fg_data com newdata
fg_surv_base <- survfit(fg_prog, newdata = perfil_base)
fg_surv_alto <- survfit(fg_prog, newdata = perfil_alto)

df_pred <- bind_rows(
  tibble(time = fg_surv_base$time / 12,
         cif  = 1 - fg_surv_base$surv,
         perfil = "M-spike = 1.0 g/dL (baixo)"),
  tibble(time = fg_surv_alto$time / 12,
         cif  = 1 - fg_surv_alto$surv,
         perfil = "M-spike = 3.0 g/dL (alto)")
)

fig5_3 <- df_pred %>%
  ggplot(aes(x = time, y = cif, color = perfil)) +
  geom_step(linewidth = 1.1) +
  scale_color_manual(values = c("#3C5488", "#E64B35"), name = NULL) +
 # scale_y_continuous(labels = scales::percent_format(accuracy = 1),
#                     limits = c(0, 0.25)) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  labs(
    title    = "Figura 5.3 — CIF de progressão ajustada: efeito do M-spike",
    subtitle = "Mulher de 70 anos. SHR do M-spike: 2.43 (Fine & Gray)",
    x        = "Anos desde diagnóstico",
    y        = "P(progressão antes do tempo t)",
    caption  = "CIF predita via Fine & Gray (finegray + coxph). mgus2, n=1384."
  )

print(fig5_3)


## ============================================================================
## DATASET 2 — pbc: CIRROSE BILIAR PRIMÁRIA — FRAGILIDADE POR ESTÁGIO
## ============================================================================
##
## Contexto:
##   Ensaio clínico randomizado (Mayo Clinic, 1974-1984) sobre D-penicilamina
##   para cirrose biliar primária. N=312 pacientes do ensaio (mais 106 fora).
##   Status: 0=censura, 1=transplante hepático, 2=morte.
##   Estágio histológico (1-4): proxy de gravidade da doença → heterogeneidade
##   clínica NÃO capturada pelas covariáveis observadas.
##
##   Referência: Fleming TR, Harrington DP (1991). Counting Processes and
##               Survival Analysis. Wiley.

pbc     <- survival::pbc 
pbc_ens <- pbc %>%
  filter(!is.na(trt)) %>%                    # apenas ensaio (n=312)
  mutate(
    morte   = as.integer(status == 2),
    log_bili = log(bili),
    stage_f  = factor(stage)
  )

cat("\n╔══════════════════════════════════════════════════════╗\n")
cat("║  DATASET: pbc — Cirrose Biliar Primária (Mayo Clinic) ║\n")
cat("╚══════════════════════════════════════════════════════╝\n\n")
cat("Status: 0=censura, 1=transplante, 2=morte\n")
print(table(pbc_ens$status, dnn = "Status"))
cat("\nEstágio histológico:\n")
print(table(pbc_ens$stage, dnn = "Estágio"))


## ----------------------------------------------------------------------------
## 5.4  Fragilidade por estágio — capturando heterogeneidade clínica
## ----------------------------------------------------------------------------
##
## TEORIA (slides 6-12):
##   µ(t | Z_estágio) = Z_estágio · µ₀(t) · exp(β₁·trt + β₂·log(bili) + β₃·albumin)
##
##   O estágio histológico define o "cluster" de fragilidade.
##   Pacientes no mesmo estágio compartilham heterogeneidade não observada
##   (microambiente tumoral, resposta imune, genética não medida).
##
##   θ > 0 indica que o estágio captura heterogeneidade ALÉM das covariáveis.

## Modelo 1: Cox padrão (sem fragilidade)
cox_pbc_std <- coxph(Surv(time, morte) ~ trt + log_bili + albumin + stage_f,
                      data = pbc_ens)

## Modelo 2: Cox com fragilidade gama por estágio
## Nota: frailty() no survival usa EM com integração analítica para gama
cox_pbc_frail <- coxph(
  Surv(time, morte) ~ trt + log_bili + albumin +
    frailty(stage, distribution = "gamma"),
  data = pbc_ens
)

cat("\n=== Comparação: Cox padrão vs. Cox com fragilidade (pbc) ===\n\n")

coef_std   <- summary(cox_pbc_std)$coefficients[1:4, c(1,2,5)]
coef_frail <- summary(cox_pbc_frail)$coefficients[1:3, c(1,2,5)]

cat("--- Modelo padrão (sem fragilidade) ---\n")
print(round(coef_std, 4))

cat("\n--- Modelo com fragilidade gama por estágio ---\n")
print(round(coef_frail, 4))

# Extrair variância da fragilidade
theta_pbc <- cox_pbc_frail$history[[1]]$history
theta_final <- theta_pbc[nrow(theta_pbc), "theta"]
theta_pbc   <- cox_pbc_frail$history[[1]]$history
theta_final <- theta_pbc[nrow(theta_pbc), "theta"]
cat("\nθ (variância da fragilidade por estágio):", round(theta_final, 4), "\n")

cat("\n--- Teste da fragilidade ---\n")
print(summary(cox_pbc_frail)$coefficients[4, ])

cat("\n--- Fragilidades estimadas por estágio ---\n")
# As fragilidades estão nos coeficientes com prefixo "gamma:"
gamma_idx_pbc <- grep("^gamma:", names(coef(cox_pbc_frail)))
frailty_vals  <- exp(coef(cox_pbc_frail)[gamma_idx_pbc])
names(frailty_vals) <- paste("Estágio", 1:4)
print(round(frailty_vals, 3))

cat("\n--- Interpretação das fragilidades por estágio ---\n")
cat("  Estágio 1: Z = 0,78  → esses pacientes têm risco 22% MENOR do que a\n")
cat("    média prevista por bili e albumina. O estágio inicial protege além\n")
cat("    do que as covariáveis observadas já capturam.\n\n")
cat("  Estágio 2: Z = 0,80  → padrão similar. Heterogeneidade residual\n")
cat("    moderada entre pacientes dentro deste estágio.\n\n")
cat("  Estágio 3: Z = 1,03  → praticamente na média. O risco clínico do\n")
cat("    estágio 3 é bem capturado por bili e albumina.\n\n")
cat("  Estágio 4: Z = 1,39  → pacientes no estágio mais avançado têm 39%\n")
cat("    a mais de risco do que as covariáveis sozinhas preveem. Há\n")
cat("    heterogeneidade não observada DENTRO do estágio 4 — pacientes com\n")
cat("    os mesmos valores de bili e albumina diferem substancialmente em\n")
cat("    suscetibilidade, possivelmente por fatores genéticos ou inflamatórios.\n\n")
cat("  θ =", round(theta_final, 4),
    "→ heterogeneidade significativa entre estágios (p = 0,023).\n")
cat("    Se θ = 0, o modelo reduziria ao Cox padrão. θ > 0 indica que\n")
cat("    o estágio captura variabilidade residual real, além do confundimento.\n\n")
cat("  Implicação: D-penicilamina tem HR ≈ 1,00 em ambos os modelos —\n")
cat("    a conclusão histórica de ineficácia do tratamento é robusta à\n")
cat("    inclusão da fragilidade por estágio.\n")

## --- Figura 5.4: Comparação de HRs --- Cox padrão vs. com fragilidade ---
hr_std   <- exp(coef(cox_pbc_std)[1:3])
hr_frail <- exp(coef(cox_pbc_frail)[1:3])
ci_std   <- exp(confint(cox_pbc_std)[1:3, ])
ci_frail <- exp(confint(cox_pbc_frail)[1:3, ])

df_hr <- bind_rows(
  tibble(
    cov    = c("Tratamento\n(D-penicilamina)", "log(Bilirrubina)", "Albumina"),
    HR     = hr_std,
    LB     = ci_std[, 1],
    UB     = ci_std[, 2],
    modelo = "Cox padrão"
  ),
  tibble(
    cov    = c("Tratamento\n(D-penicilamina)", "log(Bilirrubina)", "Albumina"),
    HR     = hr_frail,
    LB     = ci_frail[, 1],
    UB     = ci_frail[, 2],
    modelo = "Cox + fragilidade gama"
  )
) %>%
  mutate(modelo = factor(modelo, levels = c("Cox padrão", "Cox + fragilidade gama")))

fig5_4 <- df_hr %>%
  ggplot(aes(x = HR, y = cov, color = modelo, shape = modelo)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = LB, xmax = UB),
                 height = 0.15, position = position_dodge(0.5)) +
  geom_point(size = 3, position = position_dodge(0.5)) +
  scale_x_log10(breaks = c(0.1, 0.3, 0.5, 1, 2, 5, 10)) +
  scale_color_manual(values = c("#3C5488", "#E64B35"), name = "Modelo") +
  scale_shape_manual(values = c(16, 17), name = "Modelo") +
  labs(
    title    = "Figura 5.4 — Cirrose Biliar: HRs com e sem fragilidade por estágio",
    subtitle = paste("θ =", round(theta_final, 3),
                     "— heterogeneidade residual entre estágios"),
    x        = "Hazard Ratio (escala log)",
    y        = NULL,
    caption  = "pbc, Mayo Clinic (n=312). Fragilidade gama compartilhada por estágio histológico."
  )

print(fig5_4)

## DISCUSSÃO:
##   → O tratamento (D-penicilamina) não tem efeito significativo em nenhum
##     modelo. Esse resultado histórico confirmou a ineficácia do tratamento.
##   → A bilirrubina e a albumina são fortes preditores (HR >> 1 e << 1).
##   → A fragilidade por estágio captura heterogeneidade adicional: θ > 0
##     indica que pacientes no mesmo estágio ainda diferem de formas não medidas.


## ============================================================================
## DATASET 3 — veteran: CÂNCER DE PULMÃO — FRAGILIDADE POR TIPO CELULAR
## ============================================================================
##
## Contexto:
##   Ensaio do Veterans Administration (VA) sobre dois regimes de quimioterapia
##   para câncer de pulmão avançado. N=137 pacientes.
##   Tipo celular (squamous, smallcell, adeno, large) é biologicamente distinto
##   e define prognóstico. Trata-se de um cluster NATURAL com n pequeno (4 grupos).
##
##   Referência: Kalbfleisch JD, Prentice RL (1980). The Statistical Analysis
##               of Failure Time Data. Wiley.

## veteran é alias interno de "cancer" — mesmo problema do mgus2:
veteran <- survival::veteran

cat("\n╔══════════════════════════════════════════════════════╗\n")
cat("║  DATASET: veteran — Câncer de Pulmão VA (n=137)      ║\n")
cat("╚══════════════════════════════════════════════════════╝\n\n")
cat("Distribuição por tipo celular:\n")
print(table(veteran$celltype, dnn = "Tipo celular"))
cat("\nObservações: todos eventos (status=1 significa morte):\n")
print(table(veteran$status))


## ----------------------------------------------------------------------------
## 5.5  Viés de atenuação: Cox padrão vs. Cox com fragilidade
## ----------------------------------------------------------------------------
##
## O tipo celular é um fator de confusão potencial:
##   - Ele afeta a mortalidade diretamente (prognóstico biológico distinto)
##   - Ele pode estar correlacionado com outras covariáveis (idade, karno, trt)
##
## Opção 1: incluir celltype como covariável fixa (efeito fixo)
## Opção 2: tratar celltype como cluster de fragilidade (efeito aleatório)
##
## A opção 2 é mais adequada quando:
##   (a) o número de grupos é pequeno E
##   (b) nos interessa MARGINALIZAR sobre os tipos (estimativa de efeito médio)
##   (c) queremos capturar heterogeneidade residual DENTRO de cada tipo

## Modelo 1: Cox com celltype como efeito FIXO
cox_vet_fixo <- coxph(
  Surv(time, status) ~ trt + karno + age + celltype,
  data = veteran
)

## Modelo 2: Cox com celltype como fragilidade GAMA (efeito aleatório)
cox_vet_frail <- coxph(
  Surv(time, status) ~ trt + karno + age +
    frailty(celltype, distribution = "gamma"),
  data = veteran
)

## Fragilidades: coeficientes com prefixo "gamma:" precisam ser extraídos
## ANTES da tabela de HRs, pois são usados para excluir esses índices via [-gamma_idx_vet]
gamma_idx_vet <- grep("^gamma:", names(coef(cox_vet_frail)))

cat("\n=== Comparação de HRs: efeito fixo vs. fragilidade (veteran) ===\n\n")

tab_vet <- tibble(
  Covariavel      = c("Tratamento", "Score Karnofsky", "Idade"),
  HR_efeito_fixo  = round(exp(coef(cox_vet_fixo)[1:3]), 3),
  HR_fragilidade  = round(exp(coef(cox_vet_frail)[-gamma_idx_vet]), 3)
)
print(tab_vet)

## Fragilidades por tipo celular
cat("\nFragilidades estimadas (exp do efeito aleatório):\n")
z_celltype    <- exp(coef(cox_vet_frail)[gamma_idx_vet])
names(z_celltype) <- levels(veteran$celltype)
print(round(z_celltype, 3))

cat("\nVariância da fragilidade θ:\n")
hist_vet   <- cox_vet_frail$history[[1]]$history
theta_vet  <- hist_vet[nrow(hist_vet), "theta"]
cat("θ =", round(theta_vet, 4), "\n\n")

cat("--- Interpretação das fragilidades por tipo celular ---\n")
cat("  squamous: Z = 0,58  → carcinoma escamoso é 42% MENOS letal do que\n")
cat("    o previsto pelas covariáveis. Biologicamente coerente: responde\n")
cat("    melhor à quimioterapia e tem crescimento mais lento.\n\n")
cat("  large:    Z = 0,81  → abaixo da média, mas com heterogeneidade\n")
cat("    residual moderada dentro do grupo.\n\n")
cat("  smallcell: Z = 1,16 → 16% acima da média. Pequenas células são\n")
cat("    agressivas — mas parte dessa agressividade já está captada por\n")
cat("    karno (estado funcional) e idade.\n\n")
cat("  adeno:    Z = 1,45  → adenocarcinoma tem 45% mais risco residual.\n")
cat("    É o subtipo com maior heterogeneidade não observada — variação\n")
cat("    molecular intratumoral elevada que karno e idade não capturam.\n\n")

cat("--- Efeito fixo vs. fragilidade (θ =", round(theta_vet, 3), ") ---\n")
cat("  HRs de tratamento, karno e idade são muito próximos nas duas abordagens.\n")
cat("  Com k = 4 grupos pequenos, a marginalização sobre os tipos (fragilidade)\n")
cat("  vs. condicionamento (efeito fixo) produz estimativas similares.\n")
cat("  A diferença seria mais substantiva com k > 20 grupos ou θ > 0,5.\n")
cat("  A fragilidade é preferível quando há interesse na variabilidade\n")
cat("  ENTRE tipos (θ), não apenas no efeito médio dentro deles.\n")

## --- Figura 5.5: Curvas de sobrevivência por tipo celular ---
km_cell <- survfit(Surv(time, status) ~ celltype, data = veteran)
km_df_cell <- broom::tidy(km_cell) %>%
  mutate(
    tipo = gsub("celltype=", "", strata),
    time_dias = time
  )

fig5_5 <- km_df_cell %>%
  ggplot(aes(x = time_dias, y = estimate, color = tipo)) +
  geom_step(linewidth = 0.9) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = tipo),
              alpha = 0.08, color = NA) +
  scale_color_manual(
    values = c(squamous = "#3C5488", smallcell = "#E64B35",
               adeno    = "#00A087", large     = "#F39B7F"),
    name   = "Tipo celular"
  ) +
  scale_fill_manual(
    values = c(squamous = "#3C5488", smallcell = "#E64B35",
               adeno    = "#00A087", large     = "#F39B7F"),
    guide  = "none"
  ) +
  scale_x_continuous(breaks = seq(0, 1000, 200)) +
  labs(
    title    = "Figura 5.5 — Sobrevivência por tipo celular (Kaplan-Meier)",
    subtitle = paste("Heterogeneidade entre tipos: θ =", round(theta_vet, 3),
                     "— smallcell tem pior prognóstico"),
    x        = "Dias desde randomização",
    y        = "Probabilidade de sobrevivência",
    caption  = "veteran (n=137, VA trial). Bandas: IC 95% pointwise."
  )

print(fig5_5)


## ============================================================================
## SÍNTESE COMPARATIVA DOS TRÊS DATASETS
## ============================================================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  SÍNTESE — O que aprendemos com os dados reais?                     ║\n")
cat("╠══════════════════════════════════════════════════════════════════════╣\n")
cat("║                                                                     ║\n")
cat("║  mgus2 (MGUS):                                                      ║\n")
cat("║  • 1-KM superestima a progressão ao ignorar morte como competidor   ║\n")
cat("║  • Idade reduz CIF de progressão mas aumenta CS-HR → efeito         ║\n")
cat("║    indireto via competição (CS ≠ SHR quando causas interagem)       ║\n")
cat("║  • M-spike: CS ≈ SHR porque afeta apenas a progressão               ║\n")
cat("║                                                                     ║\n")
cat("║  pbc (Cirrose Biliar):                                              ║\n")
cat("║  • Fragilidade por estágio captura heterogeneidade residual (θ > 0) ║\n")
cat("║  • Mesmo com estágio como covariável, há variabilidade não explicada ║\n")
cat("║  • D-penicilamina não significativa → confirmação histórica          ║\n")
cat("║                                                                     ║\n")
cat("║  veteran (Câncer Pulmão):                                           ║\n")
cat("║  • Tipo celular como fragilidade: smallcell tem Z muito acima de 1  ║\n")
cat("║  • Marginalizar sobre tipos (fragilidade) vs. condicionar (fixo)    ║\n")
cat("║    altera a interpretação do efeito do tratamento                   ║\n")
cat("║                                                                     ║\n")
cat("╠══════════════════════════════════════════════════════════════════════╣\n")
cat("║  CONEXÃO COM A TEORIA (slides Vaupel & Yashin):                     ║\n")
cat("║  Os grupos mais frágeis (smallcell, estágio 4) morrem mais cedo,   ║\n")
cat("║  tornando os sobreviventes progressivamente mais robustos —         ║\n")
cat("║  exatamente o mecanismo do slide 4 (cáries) e slide 5 (negro/branco)║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")


## ============================================================================
## REFERÊNCIAS DOS DATASETS
## ============================================================================
##
##  mgus2:
##    Kyle RA et al. (2002). A long-term study of prognosis in monoclonal
##    gammopathy of undetermined significance. NEJM 346(8):564-569.
##    Disponível: data(mgus2, package="survival")
##
##  pbc:
##    Fleming TR, Harrington DP (1991). Counting Processes and Survival
##    Analysis. Wiley, New York.
##    D'Amico G et al. (1986). Survival and prognostic indicators in
##    compensated and decompensated cirrhosis. Hepatology 6(6):1243-1248.
##    Disponível: data(pbc, package="survival")
##
##  veteran:
##    Kalbfleisch JD, Prentice RL (1980). The Statistical Analysis of
##    Failure Time Data. Wiley, New York.
##    Disponível: data(veteran, package="survival")
##
##  Fine & Gray:
##    Fine JP, Gray RJ (1999). A proportional hazards model for the
##    subdistribution of a competing risk. JASA 94(446):496-509.
##    DOI: 10.1080/01621459.1999.10474144
##    Implementação: survival::finegray()  [Therneau TM, 2024]
##
################################################################################
