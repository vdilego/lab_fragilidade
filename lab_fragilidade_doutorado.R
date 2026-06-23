################################################################################
##                                                                            ##
##   LABORATÓRIO: HETEROGENEIDADE, FRAGILIDADE E RISCOS COMPETITIVOS          ##
##   Nível: Doutorado em Demografia — CEDEPLAR/UFMG                           ##
##   Duração estimada: 1h30                                                   ##
##   Referência principal: Vaupel, Manton & Stallard (1979); Vaupel &         ##
##   Yashin (1985); Fine & Gray (1999)                                        ##
##                                                                            ##
##   Prof. Vanessa Di Lego                                                    ##
##   Data: Junho 2026                                                         ##
##                                                                            ##
################################################################################

## ============================================================================
## ESTRUTURA DO LABORATÓRIO
## ============================================================================
##
##  PARTE 1 — Preparação e Conceitos (15 min)
##    1.1  Instalação e carregamento de pacotes
##    1.2  Simulação de população heterogênea com fragilidade gama
##    1.3  Visualizando o efeito da heterogeneidade sobre a curva observada
##
##  PARTE 2 — Modelos de Fragilidade Univariados (30 min)
##    2.1  Cox padrão (sem fragilidade): estimando o viés
##    2.2  Fragilidade gama com risco de base paramétrico (Gompertz)
##    2.3  Fragilidade gama semiparamétrica (pacote survival/frailtypack)
##    2.4  Fragilidade log-normal
##    2.5  Comparando distribuições de Z: gama vs. log-normal
##
##  PARTE 3 — Riscos Competitivos e Fine & Gray (30 min)
##    3.1  Por que o modelo de Cox simples falha com riscos competitivos
##    3.2  Incidência Cumulativa por Causa (CIF) — estimador Aalen-Johansen
##    3.3  Modelo de Fine & Gray (subdistribution hazard)
##    3.4  Fragilidade + Riscos Competitivos: o modelo de Beyersmann
##    3.5  Interpretação e armadilhas
##
##  PARTE 4 — Exercício Integrador (15 min)
##    Problema do câncer de mama (slide 44 da aula teórica)
##
## ============================================================================


## ============================================================================
## PARTE 1 — PREPARAÇÃO E CONCEITOS
## ============================================================================

## ----------------------------------------------------------------------------
## 1.1  Instalação e carregamento de pacotes
## ----------------------------------------------------------------------------

# Execute apenas na primeira vez:
# install.packages(c("survival", "frailtypack", "cmprsk", "tidycmprsk",
#                    "mstate", "tidyverse", "ggplot2", "patchwork",
#                    "flexsurv", "survminer", "ggsci"))

library(survival)
library(frailtypack)
library(cmprsk)        # Fine & Gray clássico
library(tidycmprsk)    # interface tidy para riscos competitivos
library(mstate)        # modelos multi-estado
library(tidyverse)
library(ggplot2)
library(patchwork)
library(flexsurv)      # modelos paramétricos com fragilidade
library(survminer)     # visualização de curvas de sobrevivência
library(ggsci)         # paletas de cores para publicação

# Tema ggplot para as figuras
theme_demog <- theme_bw(base_size = 13) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#2C3E50", color = NA),
        strip.text = element_text(color = "white", face = "bold"))

theme_set(theme_demog)

set.seed(2026)  # reprodutibilidade


## ----------------------------------------------------------------------------
## 1.2  Simulação: população heterogênea com fragilidade gama
## ----------------------------------------------------------------------------
##
## TEORIA (slides 6–9):
##   µ(x | Z) = Z · µ₀(x),   Z ~ Gama(1/σ², 1/σ²),   E[Z]=1, Var[Z]=σ²
##   µ₀(x) = a·exp(b·x)   [Gompertz]
##
## Vamos simular tempos de morte para N=5000 indivíduos com diferentes σ².
##

simular_populacao <- function(N = 5000, a = 0.0001, b = 0.09,
                               sigma2 = 0.5, x0 = 40) {
  # Parâmetros Gompertz: µ₀(x) = a·exp(b·x), começando na idade x0
  # Fragilidade gama: Z ~ Gama(k, k), k = 1/σ²

  k   <- 1 / sigma2
  Z   <- rgamma(N, shape = k, rate = k)   # E[Z]=1, Var[Z]=σ²

  # Tempo de morte via inversão da função de sobrevivência condicional:
  # S(x|Z) = exp(-Z · Λ₀(x)), Λ₀(x) = (a/b)(exp(bx) - 1)
  # T = (1/b) · log(1 - (b·log(U)) / (Z·a·exp(b·x0))) + x0
  # onde U ~ Uniforme(0,1)

  U <- runif(N)
  Λ0_x0 <- (a / b) * (exp(b * x0) - 1)   # risco acumulado acumulado até x0

  # Risco acumulado total que cada indivíduo deve atingir: -log(U)/Z
  H_total <- -log(U) / Z

  # Inverte Λ₀(x) = (a/b)(exp(bx)-1) para achar x
  # exp(bx) = 1 + (b/a)·H_total  →  x = log(1 + (b/a)·H_total) / b
  tempo_morte <- log(1 + (b / a) * H_total) / b

  # Censura administrativa aos 100 anos
  censurado <- as.integer(tempo_morte > 100)
  tempo_obs <- pmin(tempo_morte, 100)

  tibble(
    id      = 1:N,
    Z       = Z,
    tempo   = tempo_obs,
    evento  = 1 - censurado,   # 1 = morte observada
    sigma2  = sigma2
  )
}

# Simulamos três populações com heterogeneidade crescente
pop_homog  <- simular_populacao(sigma2 = 0.00001)  # ≈ homogênea
pop_media  <- simular_populacao(sigma2 = 0.5)
pop_alta   <- simular_populacao(sigma2 = 1.5)

pop_todas <- bind_rows(
  pop_homog  %>% mutate(grupo = "σ²≈0 (homogênea)"),
  pop_media  %>% mutate(grupo = "σ²=0.5"),
  pop_alta   %>% mutate(grupo = "σ²=1.5")
) %>% mutate(grupo = factor(grupo, levels = c("σ²≈0 (homogênea)", "σ²=0.5", "σ²=1.5")))

cat("\n=== Resumo das três populações simuladas ===\n")
pop_todas %>%
  group_by(grupo) %>%
  summarise(
    n          = n(),
    media_Z    = round(mean(Z), 4),
    var_Z      = round(var(Z), 4),
    prop_morte = round(mean(evento), 3),
    e0         = round(mean(tempo), 1)
  ) %>%
  print()


## ----------------------------------------------------------------------------
## 1.3  Visualizando o efeito da heterogeneidade
## ----------------------------------------------------------------------------

# Figura 1a: Distribuição de Z nas três populações
fig1a <- pop_todas %>%
  filter(Z < 5) %>%
  ggplot(aes(x = Z, fill = grupo, color = grupo)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_fill_jco() +
  scale_color_jco() +
  labs(title = "Distribuição da fragilidade Z",
       subtitle = "E[Z] = 1 em todas; variância difere",
       x = "Z (fragilidade individual)", y = "Densidade",
       fill = "Heterogeneidade", color = "Heterogeneidade") +
  annotate("text", x = 1.1, y = 0.1, label = "Z̄ = 1", hjust = 0, size = 3.5)

# Figura 1b: Risco observado por idade (µ̄(x)) — calculado via tabela de vida
calcular_mu_obs <- function(df, largura_intervalo = 2) {
  idades <- seq(40, 98, by = largura_intervalo)
  map_dfr(idades, function(x_inf) {
    x_sup  <- x_inf + largura_intervalo
    em_risco <- sum(df$tempo >= x_inf)
    mortes   <- sum(df$tempo >= x_inf & df$tempo < x_sup & df$evento == 1)
    mu_obs   <- mortes / (em_risco * largura_intervalo)
    tibble(idade = x_inf + largura_intervalo / 2, mu_obs = mu_obs, em_risco = em_risco)
  }) %>%
    filter(em_risco > 20)   # descarta idades com poucos indivíduos
}

mu_por_grupo <- pop_todas %>%
  group_by(grupo) %>%
  group_modify(~calcular_mu_obs(.x)) %>%
  ungroup()

# Risco de base teórico (Gompertz, a=0.0001, b=0.09)
mu_base <- tibble(
  idade  = seq(40, 98, by = 0.5),
  mu_obs = 0.0001 * exp(0.09 * idade),
  grupo  = "µ₀(x) base [Gompertz]"
)

fig1b <- mu_por_grupo %>%
  filter(mu_obs > 0) %>%
  ggplot(aes(x = idade, y = log(mu_obs), color = grupo)) +
  geom_line(linewidth = 1) +
  geom_line(data = mu_base, aes(x = idade, y = log(mu_obs)),
            color = "black", linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = 85, y = log(0.0001 * exp(0.09*85)) + 0.15,
           label = "µ₀(x) individual\n[Gompertz]", size = 3.2, hjust = 0) +
  scale_color_jco() +
  labs(title = "Risco observado µ̄(x) vs. risco individual",
       subtitle = "Heterogeneidade desacelera a curva observada",
       x = "Idade", y = "log µ̄(x)",
       color = "Heterogeneidade") +
  scale_x_continuous(breaks = seq(40, 100, 10))

fig1a + fig1b +
  plot_annotation(
    title   = "Figura 1 — Heterogeneidade e desaceleração da mortalidade",
    caption = "Baseado em Vaupel & Yashin (1985). Risco de base: µ₀(x) = 0.0001·exp(0.09x).\nN=5000 por grupo; censura aos 100 anos."
  )

## PAUSA PARA DISCUSSÃO (5 min):
##   → Por que a curva observada é mais plana que µ₀(x)?
##   → O que Vaupel & Yashin (1985) chamam de "heterogeneity's ruses"?
##   → Como isso se relaciona ao plateau de mortalidade em centenários?


## ============================================================================
## PARTE 2 — MODELOS DE FRAGILIDADE UNIVARIADOS
## ============================================================================

## Usaremos a população com σ²=0.5 como nosso "dado observado"
## (na prática, Z não é observado — só vemos tempo e evento)

dados <- pop_media %>% select(id, tempo, evento, Z)

cat("\n=== Estatísticas descritivas dos dados ===\n")
summary(dados[, c("tempo", "evento")])
cat("Proporção de mortes:", round(mean(dados$evento), 3), "\n")


## ----------------------------------------------------------------------------
## 2.1  Cox padrão (sem fragilidade) — quantificando o viés
## ----------------------------------------------------------------------------
##
## TEORIA (slide 15):
##   Ignorar fragilidade subestima os coeficientes no modelo de Cox.
##   O estimado é β_marginal ≠ β_condicional (individual).
##
## Para demonstrar o viés, adicionamos uma covariável binária X (exposição)
## que afeta Z diferencialmente.

# Criando covariável correlacionada com Z (simula variável socioeconômica)
dados <- dados %>%
  mutate(
    # Exposição: maior Z → maior probabilidade de X=1 (vulnerabilidade)
    X = rbinom(n(), 1, prob = pmin(0.9, Z / (Z + 1))),
    # Segundo grupo (controle) tem fragilidade menor em média
    Z_verdadeiro = Z
  )

# Efeito VERDADEIRO (condicional em Z): β_cond ≈ log(2) ≈ 0.69
# porque X está positivamente correlacionado com Z

# Modelo Cox SEM fragilidade
cox_sem <- coxph(Surv(tempo, evento) ~ X, data = dados)
cat("\n=== Cox sem fragilidade (β marginal — VIESADO) ===\n")
summary(cox_sem)$coefficients

# Modelo Cox COM fragilidade gama (semi-paramétrico)
cox_com_gama <- coxph(Surv(tempo, evento) ~ X + frailty(id, distribution = "gamma"),
                      data = dados)
cat("\n=== Cox com fragilidade gama (β condicional — MENOS VIESADO) ===\n")
summary(cox_com_gama)$coefficients[1, , drop = FALSE]

cat("\n--- Comparação do viés ---\n")
beta_marginal     <- coef(cox_sem)["X"]
beta_condicional  <- coef(cox_com_gama)["X"]
cat("β marginal (Cox simples):", round(beta_marginal, 4), "\n")
cat("β condicional (c/ fragilidade):", round(beta_condicional, 4), "\n")
cat("HR marginal:", round(exp(beta_marginal), 3), "\n")
cat("HR condicional:", round(exp(beta_condicional), 3), "\n")

## → O modelo sem fragilidade SUBESTIMA o HR verdadeiro.
##   Este é o "attenuation bias" por variáveis omitidas (Z).


## ----------------------------------------------------------------------------
## 2.2  Fragilidade gama com risco de base Gompertz (paramétrico)
## ----------------------------------------------------------------------------
##
## TEORIA (slide 13):
##   µ(x|Z) = Z·a·exp(bx)
##   µ̄(x) = a·exp(bx) / [1 + σ²·(a/b)·(exp(bx) - 1)]   →   plateau b/σ²
##

# Via flexsurv: especificamos Gompertz + fragilidade gama
# Nota: flexsurv usa hazard scale (shape=b, rate=a para Gompertz)
frail_gompertz <- flexsurvreg(
  Surv(tempo, evento) ~ X,
  data = dados,
  dist = "gompertz",
  mixture = FALSE
)

cat("\n=== Gompertz paramétrico (sem fragilidade explícita, para referência) ===\n")
print(frail_gompertz)

# Estimativa do plateau teórico: b/σ²
b_est    <- frail_gompertz$res["shape", "est"]
sigma2_est <- summary(cox_com_gama)$print$`Var[frailty]`
if (is.null(sigma2_est)) sigma2_est <- 0.5   # fallback

cat("\nPlateau teórico estimado (b/σ²):", round(b_est / sigma2_est, 4), "\n")
cat("Interpretação: a força de mortalidade se estabiliza em torno de",
    round(b_est / sigma2_est, 4), "por ano em idades muito avançadas.\n")


## ----------------------------------------------------------------------------
## 2.3  Fragilidade gama semiparamétrica — frailtypack
## ----------------------------------------------------------------------------
##
## frailtypack permite fragilidade compartilhada e correlacionada,
## útil para dados agrupados (famílias, clusters geográficos).
##
## Aqui usamos o modelo univariado (sem cluster real, apenas para ilustrar
## a estimação de σ²).

# Para frailtypack, precisamos de um id de cluster.
# No contexto univariado, cada indivíduo é seu próprio "cluster".
# Para um exemplo mais realista, agrupamos em 500 famílias de 10

dados_cluster <- dados %>%
  mutate(familia = rep(1:500, each = 10))

frail_gama_pack <- frailtyPenal(
  Surv(tempo, evento) ~ X + cluster(familia),
  data    = dados_cluster,
  n.knots = 8,                # nós para o risco de base não-paramétrico
  kappa   = 10000,            # penalidade de suavização
  RandDist = "Gamma"
)

cat("\n=== frailtypack: fragilidade gama compartilhada por família ===\n")
print(frail_gama_pack)

cat("\nθ (variância da fragilidade):", round(frail_gama_pack$theta, 4), "\n")
cat("IC 95%: [", round(frail_gama_pack$theta - 1.96 * frail_gama_pack$seTheta, 4),
    ",", round(frail_gama_pack$theta + 1.96 * frail_gama_pack$seTheta, 4), "]\n")

## → θ > 0 indica heterogeneidade significativa entre famílias.


## ----------------------------------------------------------------------------
## 2.4  Fragilidade log-normal
## ----------------------------------------------------------------------------
##
## TEORIA (slide 20):
##   Z = exp(W),   W ~ N(µ_Z, σ_Z²)
##   Sem forma fechada para a Laplaciana → integração numérica (EM penalizado)
##
## A log-normal tem cauda mais pesada que a gama, capturando
## heterogeneidade extrema (alguns indivíduos com Z >> 1).

frail_lognorm <- frailtyPenal(
  Surv(tempo, evento) ~ X + cluster(familia),
  data     = dados_cluster,
  n.knots  = 8,
  kappa    = 10000,
  RandDist = "LogN"
)

cat("\n=== frailtypack: fragilidade log-normal ===\n")
print(frail_lognorm)
cat("θ log-normal:", round(frail_lognorm$theta, 4), "\n")


## ----------------------------------------------------------------------------
## 2.5  Comparando distribuições: gama vs. log-normal
## ----------------------------------------------------------------------------

sigma2_gama <- frail_gama_pack$theta
sigma2_logn <- frail_lognorm$theta

z_seq <- seq(0.01, 5, length.out = 500)

dist_gama <- dgamma(z_seq, shape = 1/sigma2_gama, rate = 1/sigma2_gama)
# Log-normal: E[Z]=1 → µ = -σ²/2; σ² da log-normal ≠ θ diretamente
# θ no frailtypack para log-normal é a variância de W ~ N(0, θ)
# Então Var[Z] = (exp(theta)-1)·exp(theta) ≈ theta para theta pequeno
mu_logn   <- -sigma2_logn / 2
dist_logn <- dlnorm(z_seq, meanlog = mu_logn, sdlog = sqrt(sigma2_logn))

fig2 <- tibble(
  z      = rep(z_seq, 2),
  dens   = c(dist_gama, dist_logn),
  modelo = rep(c("Gama", "Log-normal"), each = 500)
) %>%
  ggplot(aes(x = z, y = dens, color = modelo, linetype = modelo)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "grey50") +
  scale_color_manual(values = c("#E64B35", "#4DBBD5")) +
  coord_cartesian(xlim = c(0, 4), ylim = c(0, 2)) +
  labs(
    title   = "Figura 2 — Distribuições de fragilidade estimadas",
    subtitle = "Gama vs. Log-normal: mesma média (1), caudas diferentes",
    x = "Z (fragilidade)", y = "Densidade",
    color = "Distribuição", linetype = "Distribuição",
    caption = "Parâmetros estimados dos dados simulados.\nLinhas verticais em Z=1 (média populacional)."
  )
print(fig2)

## PERGUNTA PARA DISCUSSÃO:
##   → Em que contextos a log-normal seria preferível à gama?
##   → Como o AIC/BIC pode ajudar na escolha? (veja frail_gama_pack$AIC vs frail_lognorm$AIC)
cat("\nAIC — Gama:", round(frail_gama_pack$AIC, 2),
    " | Log-normal:", round(frail_lognorm$AIC, 2), "\n")


## ============================================================================
## PARTE 3 — RISCOS COMPETITIVOS E FINE & GRAY
## ============================================================================
##
## MOTIVAÇÃO:
## Em muitos estudos de mortalidade, existem múltiplas causas de morte que
## "competem" entre si. Por exemplo:
##   Causa 1: doença cardiovascular
##   Causa 2: câncer
##   Censura: saída do estudo
##
## Tratar causas competidoras como censura produz estimativas de
## incidência cumulativa SUPER-ESTIMADAS (1 - KM não é CIF válida
## na presença de riscos competitivos).

## ----------------------------------------------------------------------------
## 3.1  Simulação com riscos competitivos
## ----------------------------------------------------------------------------

simular_riscos_competitivos <- function(N = 3000, sigma2 = 0.5) {
  # Dois riscos competitivos:
  # Causa 1 (cardiovascular): µ₁(x|Z) = Z · a1 · exp(b1·x)
  # Causa 2 (câncer):         µ₂(x|Z) = Z · a2 · exp(b2·x)
  # covariável X: aumenta risco cardiovascular, reduz câncer (plausível)

  a1 <- 0.00005; b1 <- 0.09   # cardiovascular
  a2 <- 0.00015; b2 <- 0.07   # câncer (pico mais jovem)

  k <- 1 / sigma2
  Z <- rgamma(N, shape = k, rate = k)

  X <- rbinom(N, 1, 0.5)       # exposição (ex.: tabagismo)

  # Taxa total: µ(x|Z,X) = Z·[exp(0.5·X)·µ₁(x) + exp(-0.3·X)·µ₂(x)]
  # Tempo de morte pelo risco total (inversão numérica)
  # Causa é determinada pela razão dos riscos no momento do evento

  # Simulação por inversão: tempo total
  # Λ_total(x) = Z·[exp(0.5X)·(a1/b1)(e^(b1x)-1) + exp(-0.3X)·(a2/b2)(e^(b2x)-1)]
  # Usamos o método de decomposição de causa

  U_total <- runif(N)
  # Tentativa iterativa (Newton simples) para resolver Λ_total(x) = -log(U)
  tempo_morte <- numeric(N)
  x_iter <- rep(50, N)   # início em 50 anos

  Lambda_total <- function(x, z, xi) {
    z * (exp(0.5 * xi) * (a1/b1) * (exp(b1*x) - 1) +
         exp(-0.3 * xi) * (a2/b2) * (exp(b2*x) - 1))
  }
  lambda_total <- function(x, z, xi) {
    z * (exp(0.5 * xi) * a1 * exp(b1*x) +
         exp(-0.3 * xi) * a2 * exp(b2*x))
  }
  alvo <- -log(U_total)

  for (iter in 1:50) {
    f_val <- Lambda_total(x_iter, Z, X) - alvo
    f_der <- lambda_total(x_iter, Z, X)
    x_iter <- pmax(0, x_iter - f_val / f_der)
  }
  tempo_morte <- x_iter

  # Determinar causa: probabilidade de ser causa 1 no momento do evento
  mu1_t <- exp(0.5 * X) * a1 * exp(b1 * tempo_morte)
  mu2_t <- exp(-0.3 * X) * a2 * exp(b2 * tempo_morte)
  prob_causa1 <- mu1_t / (mu1_t + mu2_t)
  causa <- ifelse(runif(N) < prob_causa1, 1, 2)

  # Censura administrativa aos 90 anos
  censurado    <- tempo_morte > 90
  tempo_obs    <- pmin(tempo_morte, 90)
  status       <- ifelse(censurado, 0, causa)   # 0=censura, 1=cardio, 2=câncer

  tibble(
    id      = 1:N,
    Z       = Z,
    tempo   = tempo_obs,
    status  = status,                           # 0, 1 ou 2
    evento1 = as.integer(status == 1),          # indicador causa 1
    evento2 = as.integer(status == 2),          # indicador causa 2
    X       = X,
    sigma2  = sigma2
  )
}

dados_cr <- simular_riscos_competitivos(N = 3000, sigma2 = 0.5)

cat("\n=== Distribuição dos desfechos (riscos competitivos) ===\n")
table(dados_cr$status, dnn = "Status (0=censura, 1=cardio, 2=câncer)")


## ----------------------------------------------------------------------------
## 3.2  Incidência Cumulativa por Causa — CIF (Aalen-Johansen)
## ----------------------------------------------------------------------------
##
## F_k(t) = P(T ≤ t, causa = k) = ∫₀ᵗ S(u⁻) · hk(u) du
##
## Esta é a estimativa CORRETA de incidência com riscos competitivos.
## NÃO é 1 - Kaplan-Meier (que trata outras causas como censura aleatória
## e superestima a incidência de cada causa).

# Estimativa KM (ERRADA para incidência com riscos competitivos)
km_causa1_errado <- survfit(Surv(tempo, evento1) ~ X, data = dados_cr)

# CIF correta via cmprsk
cif_obj <- with(dados_cr, cuminc(ftime = tempo, fstatus = status,
                                  group = X, cencode = 0))

# Extraindo CIF para plotar
extrair_cif <- function(obj, grupo, causa) {
  nome <- paste(grupo, causa)
  tibble(
    tempo    = obj[[nome]]$time,
    cif      = obj[[nome]]$est,
    grupo    = as.character(grupo),
    causa    = paste("Causa", causa)
  )
}

df_cif <- bind_rows(
  extrair_cif(cif_obj, 0, 1), extrair_cif(cif_obj, 0, 2),
  extrair_cif(cif_obj, 1, 1), extrair_cif(cif_obj, 1, 2)
) %>%
  mutate(grupo = ifelse(grupo == "0", "X=0 (não exposto)", "X=1 (exposto)"))

# Comparação KM (errado) vs. CIF (correto) para causa 1
km_df <- broom::tidy(km_causa1_errado) %>%
  mutate(
    incidencia_km = 1 - estimate,
    grupo = ifelse(strata == "X=0", "X=0 (não exposto)", "X=1 (exposto)")
  )

fig3a <- ggplot() +
  geom_step(data = km_df,
            aes(x = time, y = incidencia_km, color = grupo, linetype = "1-KM (ERRADO)"),
            linewidth = 0.9) +
  geom_step(data = df_cif %>% filter(causa == "Causa 1"),
            aes(x = tempo, y = cif, color = grupo, linetype = "CIF (CORRETO)"),
            linewidth = 0.9) +
  scale_color_manual(values = c("#E64B35", "#3C5488")) +
  scale_linetype_manual(values = c("dashed", "solid"),
                        name = "Estimador") +
  labs(
    title   = "Figura 3a — 1-KM vs. CIF para Causa 1 (cardiovascular)",
    subtitle = "1-KM superestima a incidência quando há riscos competitivos",
    x = "Idade", y = "Probabilidade de incidência",
    color = "Grupo"
  )

fig3b <- df_cif %>%
  ggplot(aes(x = tempo, y = cif, color = grupo, linetype = causa)) +
  geom_step(linewidth = 0.9) +
  scale_color_manual(values = c("#E64B35", "#3C5488")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(
    title   = "Figura 3b — CIF por causa e grupo",
    subtitle = "Incidências cumulativas somam S(t) + F₁(t) + F₂(t) = 1",
    x = "Idade", y = "CIF",
    color = "Grupo", linetype = "Causa"
  )

(fig3a / fig3b) +
  plot_annotation(caption = "cmprsk::cuminc(). Dados simulados: N=3000, σ²=0.5.")

## PONTO CRÍTICO:
##   A soma das duas CIFs + S(t) deve ser igual a 1 para cada grupo.
##   Isso é impossível se usarmos 1-KM para cada causa separadamente!

df_cif %>%
  group_by(grupo) %>%
  filter(tempo == max(tempo)) %>%
  group_by(grupo, causa) %>%
  summarise(max_cif = round(max(cif), 3)) %>%
  pivot_wider(names_from = causa, values_from = max_cif) %>%
  mutate(soma_CIFs = `Causa 1` + `Causa 2`) %>%
  print()


## ----------------------------------------------------------------------------
## 3.3  Modelo de Fine & Gray — subdistribution hazard
## ----------------------------------------------------------------------------
##
## REFERÊNCIA: Fine JP, Gray RJ (1999). A proportional hazards model for
##             the subdistribution of a competing risk. JASA 94(446):496-509.
##
## IDEIA CENTRAL:
##   Ao invés de modelar o hazard de causa específica h_k(t) (abordagem de
##   causa-específica), Fine & Gray modelam o "subdistribution hazard":
##
##   h̃_k(t) = -d/dt [log(1 - F_k(t))]
##
##   Interpretação: taxa instantânea de evento k para indivíduos que AINDA
##   NÃO tiveram o evento k — INCLUINDO os que já morreram por outra causa
##   (eles permanecem em "risco" para k com weight → 0).
##
##   Propriedade fundamental: exp(β̃_k) descreve diretamente o efeito de
##   uma covariável sobre F_k(t). Portanto é o modelo correto quando o
##   interesse é na incidência cumulativa, não no hazard condicional.
##
##   CONTRASTE com causa-específica:
##   • Causa-específica (Cox por causa): "dada sobrevivência, qual é o risco?"
##     → útil para mecanismo etiológico
##   • Fine & Gray: "qual é o impacto na incidência acumulada observada?"
##     → útil para prognóstico e comunicação de risco ao paciente
##

# Fine & Gray clássico via cmprsk
fg_causa1 <- crr(
  ftime   = dados_cr$tempo,
  fstatus = dados_cr$status,
  cov1    = model.matrix(~ X, data = dados_cr)[, -1, drop = FALSE],
  failcode = 1,
  cencode  = 0
)

fg_causa2 <- crr(
  ftime   = dados_cr$tempo,
  fstatus = dados_cr$status,
  cov1    = model.matrix(~ X, data = dados_cr)[, -1, drop = FALSE],
  failcode = 2,
  cencode  = 0
)

cat("\n=== Fine & Gray — Subdistribution Hazard Ratio ===\n")
cat("\n--- Causa 1 (cardiovascular) ---\n")
summary(fg_causa1)

cat("\n--- Causa 2 (câncer) ---\n")
summary(fg_causa2)

# Comparação com causa-específica (Cox por causa)
cox_causa1 <- coxph(Surv(tempo, evento1) ~ X, data = dados_cr)
cox_causa2 <- coxph(Surv(tempo, evento2) ~ X, data = dados_cr)

cat("\n=== Comparação: Cause-Specific HR vs. Subdistribution HR ===\n")
tab_comp <- tibble(
  Abordagem = c("Causa-específica (Cox)", "Fine & Gray (subdist.)"),
  HR_causa1 = c(
    round(exp(coef(cox_causa1)["X"]), 3),
    round(exp(fg_causa1$coef["X"]), 3)
  ),
  HR_causa2 = c(
    round(exp(coef(cox_causa2)["X"]), 3),
    round(exp(fg_causa2$coef["X"]), 3)
  )
)
print(tab_comp)

cat("\nInterpretação:\n")
cat("• CS-HR causa 1: efeito de X sobre o risco de morte cardiovascular\n")
cat("  DADO QUE o indivíduo está vivo (não morreu por outra causa)\n")
cat("• SHR causa 1: efeito de X sobre a incidência cumulativa cardiovascular\n")
cat("  (inclui o efeito indireto via competição com causa 2)\n")

## NOTA IMPORTANTE: quando as duas causas competem, o SHR e o CS-HR divergem.
## Se X aumenta risco de causa 1 E reduz risco de causa 2,
## o SHR de causa 1 é MAIOR que o CS-HR (porque a exposição mantém
## mais pessoas em risco de causa 1 ao reduzir morte por causa 2).


## Visualização: CIF ajustada por Fine & Gray
# Usando tidycmprsk para interface mais moderna
library(tidycmprsk)

dados_cr_factor <- dados_cr %>%
  mutate(status_f = factor(status, levels = c(0, 1, 2),
                           labels = c("censura", "cardio", "cancer")),
         X_f = factor(X, labels = c("Não exposto", "Exposto")))

fg_tidy_c1 <- tidycmprsk::crr(
  Surv(tempo, status_f) ~ X_f,
  data     = dados_cr_factor,
  failcode = "cardio"
)

cat("\n=== Fine & Gray via tidycmprsk (causa 1) ===\n")
tbl_regression(fg_tidy_c1, exponentiate = TRUE)


## ----------------------------------------------------------------------------
## 3.4  Fragilidade + Riscos Competitivos
## ----------------------------------------------------------------------------
##
## Extensão de Beyersmann et al. (2012) e Rondeau et al. (2015):
## Combina fragilidade compartilhada com subdistribution hazard para
## dados agrupados com riscos competitivos.
##
## O modelo: h̃_k(t|Z) = Z · h̃_k0(t) · exp(X'β_k)
##
## Implementação via frailtypack::multivPenal

dados_cr_cluster <- dados_cr %>%
  mutate(familia = rep(1:300, each = 10))

cat("\n=== Fragilidade + Riscos Competitivos (frailtypack::multivPenal) ===\n")
tryCatch({
  frail_cr <- multivPenal(
    Surv(tempo, status) ~ X + cluster(familia),
    data       = dados_cr_cluster,
    n.knots    = 6,
    kappa      = c(1000, 1000),
    RandDist   = "Gamma",
    competing.risks = TRUE
  )
  print(frail_cr)
  cat("θ (fragilidade compartilhada nos riscos competitivos):",
      round(frail_cr$theta, 4), "\n")
}, error = function(e) {
  cat("Nota: multivPenal requer versão >= 3.3 do frailtypack.\n")
  cat("Alternativa: usar cause-specific frailty Cox separadamente.\n")
  cat("  cox_cr1 <- coxph(Surv(tempo, evento1) ~ X + frailty(familia), data=dados_cr_cluster)\n")
  cat("  cox_cr2 <- coxph(Surv(tempo, evento2) ~ X + frailty(familia), data=dados_cr_cluster)\n")
})


## ----------------------------------------------------------------------------
## 3.5  Interpretação e armadilhas — Resumo comparativo
## ----------------------------------------------------------------------------

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║        GUIA DE ESCOLHA: QUAL MODELO USAR?                           ║\n")
cat("╠══════════════════════════════════════════════════════════════════════╣\n")
cat("║ Pergunta de pesquisa        │ Modelo adequado                       ║\n")
cat("╠══════════════════════════════════════════════════════════════════════╣\n")
cat("║ Mecanismo etiológico        │ Cox causa-específico                   ║\n")
cat("║ (dado vivo, qual o risco?)  │ (uma regressão por causa)             ║\n")
cat("╠══════════════════════════════════════════════════════════════════════╣\n")
cat("║ Prognóstico / incidência    │ Fine & Gray (SHR)                     ║\n")
cat("║ observada (o que acontece   │ Modela F_k(t) diretamente             ║\n")
cat("║ na população?)              │                                       ║\n")
cat("╠══════════════════════════════════════════════════════════════════════╣\n")
cat("║ Dados agrupados / famílias  │ Fragilidade compartilhada             ║\n")
cat("║ (dependência intracluster)  │ (frailtypack, coxme)                  ║\n")
cat("╠══════════════════════════════════════════════════════════════════════╣\n")
cat("║ Agrupados + riscos compet.  │ multivPenal / cause-specific frailty  ║\n")
cat("╠══════════════════════════════════════════════════════════════════════╣\n")
cat("║ Plateau / centenários       │ Fragilidade + Gompertz paramétrico    ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n")


## ============================================================================
## PARTE 4 — EXERCÍCIO INTEGRADOR
## ============================================================================
##
## Baseado no EXERCÍCIO DO SLIDE 44 da aula teórica:
##
## "Uma pesquisadora estuda mortalidade por câncer de mama (1970-2010).
##  Após ajustar por estágio e tratamento, a mortalidade de baixa renda
##  começa mais alta mas CRUZA a de alta renda por volta dos 75 anos."
##
## Simulamos esse cenário e pedimos que os/as doutorandos/as respondam
## as 5 perguntas do slide usando as ferramentas do laboratório.

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║  EXERCÍCIO INTEGRADOR — Câncer de mama com riscos competitivos      ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

## Simulação do cenário do exercício
## Baixa renda (renda=0): maior fragilidade média (mais heterogênea)
## Alta renda  (renda=1): menor fragilidade média (mais homogênea)

set.seed(42)
N_ex <- 2000

# Grupo de baixa renda: σ²=1.2 (alta heterogeneidade → seleção intensa)
baixa_renda <- simular_populacao(N = N_ex, a = 0.0002, b = 0.08,
                                  sigma2 = 1.2, x0 = 30) %>%
  mutate(renda = 0, grupo = "Baixa renda")

# Grupo de alta renda: σ²=0.2 (baixa heterogeneidade → seleção fraca)
alta_renda  <- simular_populacao(N = N_ex, a = 0.00008, b = 0.08,
                                  sigma2 = 0.2, x0 = 30) %>%
  mutate(renda = 1, grupo = "Alta renda")

dados_ex <- bind_rows(baixa_renda, alta_renda)

## --- Pergunta 1: Mostrar o cruzamento das curvas ---

cat("PERGUNTA 1: Quais explicações são possíveis para o cruzamento?\n\n")

# Calcular µ observado por grupo
mu_ex <- dados_ex %>%
  group_by(grupo) %>%
  group_modify(~calcular_mu_obs(.x, largura_intervalo = 3)) %>%
  ungroup()

fig_ex1 <- mu_ex %>%
  filter(mu_obs > 0, em_risco > 30) %>%
  ggplot(aes(x = idade, y = log(mu_obs), color = grupo)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_manual(values = c("#E64B35", "#3C5488"),
                     name = "Grupo socioeconômico") +
  labs(
    title   = "EXERCÍCIO — Cruzamento das curvas de mortalidade",
    subtitle = paste("Baixa renda: σ²=1.2 (alta heterog.) | Alta renda: σ²=0.2 (baixa heterog.)\n",
                     "Cruzamento em torno dos 70-75 anos é artefato de seleção diferencial"),
    x = "Idade", y = "log µ̄(x)",
    caption = "N=2000 por grupo. O cruzamento é produzido puramente pela diferença em σ², não por mudança biológica."
  ) +
  scale_x_continuous(breaks = seq(30, 100, 10))

print(fig_ex1)

## --- Pergunta 2: Escrever o modelo de fragilidade ---
cat("\nPERGUNTA 2: Modelo de fragilidade univariado\n")
cat("  µ(x|Z, renda) = Z · µ₀(x) · exp(β · renda)\n")
cat("  Z ~ Gama(1/σ², 1/σ²)\n")
cat("  µ₀(x) = a·exp(bx)  [Gompertz]\n")
cat("  Parâmetros: a, b, β, σ²\n")

# Estimação
dados_ex_cox <- dados_ex %>%
  mutate(familia_ex = rep(1:400, each = 10))

cox_ex_simples <- coxph(Surv(tempo, evento) ~ renda, data = dados_ex)
cox_ex_frail   <- coxph(Surv(tempo, evento) ~ renda +
                         frailty(familia_ex, distribution = "gamma"),
                        data = dados_ex)

cat("\nHR(renda) sem fragilidade:", round(exp(coef(cox_ex_simples)["renda"]), 3), "\n")
cat("HR(renda) com fragilidade:", round(exp(coef(cox_ex_frail)["renda"]), 3), "\n")

## --- Pergunta 3: Identificabilidade ---
cat("\nPERGUNTA 3: Identificabilidade\n")
cat("  Sem covariáveis: NÃO identificável (Elbers & Ridder, 1982)\n")
cat("  COM covariável 'renda': IDENTIFICÁVEL (E[Z] < ∞ garantido pela gama)\n")
cat("  Seria necessário: especificar µ₀(x) parametricamente OU\n")
cat("  ter covariável com efeito conhecido (instrumento)\n")

## --- Pergunta 4: Viés no Cox ---
cat("\nPERGUNTA 4: Viés de covariáveis ao ignorar fragilidade\n")
cat("  HR marginal (sem fragilidade):", round(exp(coef(cox_ex_simples)["renda"]), 3), "\n")
cat("  HR condicional (com fragilidade):", round(exp(coef(cox_ex_frail)["renda"]), 3), "\n")
cat("  → O modelo sem fragilidade SUBESTIMA o efeito protetor da alta renda.\n")
cat("  → Isso ocorre porque baixa renda tem maior seleção de frágeis:\n")
cat("    os sobreviventes de baixa renda são progressivamente mais robustos.\n")

## --- Pergunta 5: Dados para fragilidade correlacionada ---
cat("\nPERGUNTA 5: Fragilidade correlacionada\n")
cat("  Dado adicional: histórico familiar de câncer de mama (pares de irmãs/mães-filhas)\n")
cat("  Modelo: (Z₁, Z₂) com Cor(Z₁, Z₂) = ρ — não compartilham Z idêntico\n")
cat("  Nova pergunta: Quanto da variação em suscetibilidade é genética vs. ambiental?\n")
cat("  ρMZ = correlação em gêmeas monozigóticas\n")
cat("  ρDZ = correlação em gêmeas dizigóticas\n")
cat("  Herdabilidade ≈ 2(ρMZ - ρDZ)  [fórmula de Falconer]\n")


## ============================================================================
## FIGURA FINAL: SÍNTESE DO LABORATÓRIO
## ============================================================================

cat("\n=== Gerando figura síntese ===\n")

# Evolução da fragilidade média residual Z̄(x) com a idade
# Z̄(x) = 1 / [1 + σ²·Λ₀(x)]   [fórmula do slide 9]

calcular_Zbar <- function(x, sigma2, a = 0.0001, b = 0.09) {
  Lambda0 <- (a/b) * (exp(b*x) - 1)
  1 / (1 + sigma2 * Lambda0)
}

fig_sintese <- expand_grid(
  idade  = seq(40, 100, by = 1),
  sigma2 = c(0.1, 0.5, 1.0, 1.5)
) %>%
  mutate(
    Zbar  = calcular_Zbar(idade, sigma2),
    label = paste("σ² =", sigma2)
  ) %>%
  ggplot(aes(x = idade, y = Zbar, color = factor(sigma2))) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey50") +
  scale_color_jco(name = "σ² (heterog.)") +
  annotate("text", x = 42, y = 1.01, label = "Z̄ inicial = 1", size = 3.2, hjust = 0) +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(
    title   = "Figura Síntese — Seleção de sobreviventes: queda de Z̄(x) com a idade",
    subtitle = "Quanto maior σ², mais rápido Z̄(x) cai: os frágeis morrem cedo, os robustos sobrevivem",
    x = "Idade",
    y = "Fragilidade média dos sobreviventes Z̄(x)",
    caption = paste("Fórmula: Z̄(x) = 1/[1 + σ²·Λ₀(x)], com µ₀(x) = 0.0001·exp(0.09x).",
                    "\nVaupel, Manton & Stallard (1979); Vaupel & Yashin (1985).")
  )

print(fig_sintese)


## ============================================================================
## REFERÊNCIAS
## ============================================================================
##
##  Vaupel JW, Manton KG, Stallard E (1979). The impact of heterogeneity in
##    individual frailty on the dynamics of mortality. Demography 16(3):439-454.
##
##  Vaupel JW, Yashin AI (1985). Heterogeneity's ruses: some surprising
##    effects of selection on population dynamics. Am. Statistician 39(3):176-185.
##
##  Vaupel JW, Yashin AI (1987). Repeated resuscitation: how lifesaving
##    alters life tables. Demography 24(1):123-135.
##
##  Vaupel JW, Yashin AI (1987). Targeting lifesaving: demographic linkages
##    between population structure and life expectancy. Eur. J. Population 2(3):335-360.
##
##  Vaupel JW, Yashin AI, Manton KG (1988). Debilitation's aftermath:
##    stochastic process models of mortality. Math. Population Studies 1(1):21-48.
##
##  Elbers C, Ridder G (1982). True and spurious duration dependence: the
##    identifiability of the proportional hazard model. Rev. Econ. Studies 49(3):403-409.
##
##  Fine JP, Gray RJ (1999). A proportional hazards model for the subdistribution
##    of a competing risk. JASA 94(446):496-509.
##
##  Wienke A (2003). Frailty models. MPIDR Working Paper WP 2003-032.
##
##  Duchateau L, Janssen P (2008). The Frailty Model. Springer.
##
##  Beyersmann J, Allignol A, Schumacher M (2012). Competing Risks and
##    Multistate Models with R. Springer.
##
##  Rondeau V, Gonzalez JR, Mazroui Y, Mauguen A, Diakite A, Laurent A,
##    Dumerc M, Krol A, Calvo E (2024). frailtypack: Shared Frailty Models,
##    Nested Frailty Models, Frailty Models for Joint Outcomes. R package.
##
## ============================================================================

cat("\n\nLaboratório concluído.\n")
cat("Para dúvidas ou extensões, consulte: vdilego.github.io\n")
