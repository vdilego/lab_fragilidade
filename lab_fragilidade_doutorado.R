################################################################################
##                                                                            ##
##   LABORATÓRIO: HETEROGENEIDADE, FRAGILIDADE E RISCOS COMPETITIVOS          ##
##   Referência principal: Vaupel, Manton & Stallard (1979); Vaupel &         ##
##   Yashin (1985); Fine & Gray (1999)                                        ##
##                                                                            ##
##   Prof. Vanessa Di Lego                                                    ##
##   Data: Junho 2026                                                         ##
##                                                                            ##
################################################################################


## ============================================================================
## PARTE 1 — PREPARAÇÃO E CONCEITOS
## ============================================================================

## ----------------------------------------------------------------------------
## 1.1  Instalação e carregamento de pacotes
## ----------------------------------------------------------------------------

# Execute apenas na primeira vez:
# install.packages(c("survival", "cmprsk", "tidyverse",
#                    "ggplot2", "patchwork", "ggsci", "scales", "broom"))
#
# Pacotes opcionais (instalação separada, não necessários para este script):
#   frailtypack  — fragilidade semiparamétrica avançada
#   flexsurv     — modelos paramétricos com fragilidade gama
#   survminer    — visualizações Kaplan-Meier prontas para publicação
#   tidycmprsk   — interface tidy para riscos competitivos

library(survival)   # coxph, frailty, finegray, survfit — núcleo do laboratório
library(cmprsk)     # Fine & Gray via crr() — alternativa ao finegray()
library(tidyverse)  # dplyr, ggplot2, purrr, tibble
library(ggplot2)
library(patchwork)  # combinar múltiplos gráficos
library(ggsci)      # paletas de cores para publicação (jco, nejm)
library(broom)      # tidy() para objetos survfit e coxph

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

## Parâmetros globais — mesmos valores dos defaults da função acima.
## Definidos aqui (após as simulações) para que as seções 2.2 e seguintes
## possam usar a, b, sigma2 e x0 diretamente, sem redefini-los.
## IMPORTANTE: estes valores devem ser consistentes com a chamada
##   simular_populacao(sigma2 = 0.5) usada em `pop_media` acima.
a      <- 0.0001   # nível base do risco Gompertz: µ₀(x) = a·exp(b·x)
b      <- 0.09     # ritmo de aumento do risco com a idade
sigma2 <- 0.5      # variância da fragilidade (população de referência)
x0     <- 40       # idade de início da observação

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

## RESULTADO CONTRAINTUITIVO — LEIA ANTES DE CONTINUAR:
##
##   A tabela acima mostra que e₀ SOBE com σ² (ex.: ~69 → ~72 → ~79 anos).
##   Isso parece contraditório: indivíduos com Z alto têm vida mais curta,
##   então por que a média da população aumenta com a heterogeneidade?
##
##   Resposta: DESIGUALDADE DE JENSEN.
##
##   e₀(Z) é uma função CONVEXA de Z:
##     e₀(Z=0.25) ≈ 85 anos  (+14 acima de e₀(1))
##     e₀(Z=1.00) ≈ 71 anos  [referência]
##     e₀(Z=1.75) ≈ 66 anos  (−5 abaixo de e₀(1))
##
##   O ganho dos robustos (+14) supera a perda dos frágeis (−5).
##   Formalmente: E[e₀(Z)] ≥ e₀(E[Z]) = e₀(1) para qualquer f convexa.
##   Maior σ² → distribuição mais espalhada → maior ganho Jensen → maior e₀.
##
##   Isso NÃO contradiz o plateau: µ̄(x) desacelera (plateau mais baixo com
##   maior σ²) ao mesmo tempo que e₀ sobe. São dois fenômenos distintos:
##   o plateau descreve a COMPOSIÇÃO dos sobreviventes; e₀ descreve a
##   MÉDIA da população original.

cat("\n=== Verificação da convexidade de e₀(Z) ===\n")
e0_dado_Z <- function(Z_val, a=0.0001, b=0.09, x0=40) {
  # e₀(Z) = x0 + integral de S(x|Z) dx
  # S(x|Z) = exp(-Z * (Lambda0(x) - Lambda0(x0)))
  Lambda0 <- function(x) (a/b)*(exp(b*x) - 1)
  L0x0    <- Lambda0(x0)
  x0 + integrate(function(x) exp(-Z_val*(Lambda0(x) - L0x0)), x0, x0+200)$value
}
tab_e0Z <- data.frame(
  Z       = c(0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 3.00),
  e0_anos = sapply(c(0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 3.00), e0_dado_Z)
)
tab_e0Z$var_vs_Z1 <- round(tab_e0Z$e0_anos - tab_e0Z$e0_anos[tab_e0Z$Z == 1], 1)
tab_e0Z$e0_anos   <- round(tab_e0Z$e0_anos, 1)
print(tab_e0Z)

cat("\nAssimetria Jensen: ganho de Z=0.25 (+", tab_e0Z$var_vs_Z1[1],
    " anos) > perda de Z=1.75 (",
    round(e0_dado_Z(1.75)-e0_dado_Z(1.0), 1), " anos)\n")
cat("→ A média ponderada E[e₀(Z)] sempre supera e₀(E[Z]) = e₀(1).\n")


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
##
##   Sobre µ̄(x) e o plateau:
##   → Por que a curva observada é mais plana que µ₀(x)?
##   → O que Vaupel & Yashin (1985) chamam de "heterogeneity's ruses"?
##   → Como isso se relaciona ao plateau de mortalidade em centenários?
##
##   Sobre e₀ e a desigualdade de Jensen:
##   → A tabela acima mostra e₀ SUBINDO com σ². Por que isso não contradiz
##     o fato de que indivíduos com Z alto morrem mais cedo?
##   → O que é uma função convexa? Por que e₀(Z) é convexa?
##   → Como reconciliar "plateau mais baixo com σ² maior" (µ̄(x) desacelera)
##     com "e₀ maior com σ² maior" (Jensen)?
##   → Se σ² → ∞, o que acontece com e₀? E com o plateau b/σ²?


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
# e agrupando em clusters de família para a demonstração de fragilidade.
#
# NOTA TÉCNICA: frailty(id) com um cluster por indivíduo (N=5000 clusters)
#   produz θ → 0 por overfitting e pode não convergir. O correto para
#   demonstrar o viés é usar clusters com múltiplos indivíduos.
#   Aqui usamos 500 famílias de 10 indivíduos cada.

dados <- dados %>%
  mutate(
    X            = rbinom(n(), 1, prob = pmin(0.9, Z / (Z + 1))),
    Z_verdadeiro = Z,
    familia      = as.factor(rep(1:500, each = 10))
  )

# Efeito VERDADEIRO (condicional em Z): β_cond ≈ log(2) ≈ 0.69
# porque X está positivamente correlacionado com Z

# Modelo Cox SEM fragilidade — estima β_marginal (populacional)
cox_sem <- coxph(Surv(tempo, evento) ~ X, data = dados)
cat("\n=== Cox sem fragilidade (β marginal — VIESADO) ===\n")
print(round(summary(cox_sem)$coefficients, 4))

# Modelo Cox COM fragilidade gama — estima β_condicional (individual)
# Usa clusters de família (500 × 10) para identificabilidade adequada.
cox_com_gama <- coxph(
  Surv(tempo, evento) ~ X + frailty(familia, distribution = "gamma"),
  data = dados
)
cat("\n=== Cox com fragilidade gama (β condicional — MENOS VIESADO) ===\n")
print(round(summary(cox_com_gama)$coefficients[1, , drop = FALSE], 4))

cat("\n--- Comparação do viés ---\n")
beta_marginal    <- coef(cox_sem)["X"]
beta_condicional <- coef(cox_com_gama)["X"]
cat("β marginal (Cox simples):       ", round(beta_marginal, 4), "\n")
cat("β condicional (c/ fragilidade): ", round(beta_condicional, 4), "\n")
cat("HR marginal:    ", round(exp(beta_marginal), 3), "\n")
cat("HR condicional: ", round(exp(beta_condicional), 3), "\n")
cat("Viés (subestimativa):", round((1 - exp(beta_marginal)/exp(beta_condicional))*100, 1),
    "%\n")

## → O modelo sem fragilidade SUBESTIMA o HR verdadeiro.
##   Este é o "attenuation bias" por variáveis omitidas (Z).
##   O HR marginal é sempre menor que o HR condicional quando X e Z
##   estão positivamente correlacionados.


## ----------------------------------------------------------------------------
## 2.2  Plateau de Gompertz + fragilidade gama: verificação analítica e numérica
## ----------------------------------------------------------------------------
##
## TEORIA (slide 13):
##   µ(x|Z) = Z · a · exp(bx)
##   µ̄(x)  = a · exp(bx) / [1 + σ² · (a/b) · (exp(bx) − 1)]
##   Limite para x → ∞:   µ̄(x) → b/σ²   [plateau]
##
## NOTA: {flexsurv} não é necessário aqui. O plateau b/σ² decorre da
##   matemática do modelo e pode ser:
##   (a) calculado diretamente com os parâmetros conhecidos da simulação
##   (b) estimado empiricamente por regressão de log(µ̄(x)) sobre idade
##   (c) verificado numericamente comparando µ̄(x) observada com b/σ²
##
## As três abordagens são mostradas abaixo.

## --- (a) Plateau analítico com parâmetros da simulação ---
# Os parâmetros a, b e sigma2 foram definidos na simulação acima.
plateau_teorico <- b / sigma2   # b=0.09, sigma2=0.5 → plateau=0.18

cat("\n=== Plateau de Gompertz + fragilidade gama ===\n")
cat("Parâmetros da simulação: a =", a, " b =", b, " σ² =", sigma2, "\n")
cat("Plateau teórico b/σ² =", round(plateau_teorico, 4),
    "por ano (µ̄(x) converge para este valor em x→∞)\n")

## --- (b) Estimativa de b por regressão de log(µ̄(x)) sobre idade ---
## Nas idades jovens (antes da seleção distorcer a curva observada),
## log(µ̄(x)) ≈ log(a) + bx porque Z̄(x) ≈ 1.
## Usamos apenas idades 40–65 para a regressão.
calcular_mu_obs_2 <- function(df, largura = 5) {
  idades <- seq(40, 95, by = largura)
  res <- lapply(idades, function(xi) {
    n <- sum(df$tempo >= xi)
    d <- sum(df$tempo >= xi & df$tempo < xi + largura & df$evento == 1)
    if (n < 50 || d == 0) return(NULL)
    data.frame(idade = xi + largura/2, mu_obs = d / (n * largura), n_risco = n)
  })
  do.call(rbind, Filter(Negate(is.null), res))
}

mu_tabela <- calcular_mu_obs_2(dados)

# Regressão apenas nas idades jovens (antes da seleção)
mu_jovem <- subset(mu_tabela, idade < 66)
lm_gomp  <- lm(log(mu_obs) ~ idade, data = mu_jovem)
b_est    <- coef(lm_gomp)["idade"]
a_est    <- exp(coef(lm_gomp)["(Intercept)"])

cat("\nGompertz estimado por regressão de log(µ̄(x)) sobre idade (idades 40–65):\n")
cat("  a_est =", round(a_est, 6), "  b_est =", round(b_est, 4), "\n")
cat("  Valores verdadeiros: a = 0.0001  b = 0.09\n")
cat("  Plateau estimado b_est/σ²:", round(b_est / sigma2, 4), "\n")

## --- (c) Verificação numérica: µ̄(x) observada se aproxima de b/σ²? ---
cat("\n--- Verificação numérica: µ̄(x) nas idades avançadas ---\n")
cat("(compara com plateau teórico =", round(plateau_teorico, 5), ")\n\n")
print(round(tail(mu_tabela, 5), 5))

cat("\nInterpretação:\n")
cat("  µ̄(x) observada nas faixas 90–100 anos se aproxima de", round(plateau_teorico, 3), "\n")
cat("  mas nunca o atinge completamente com dados finitos.\n")
cat("  Maior σ² → plateau MAIS BAIXO (b/σ²↓); µ̄(x) desacelera mais cedo.\n")
cat("  A curva individual µ₀(x) = a·exp(bx) nunca para de crescer —\n")
cat("  o plateau é propriedade da POPULAÇÃO, não do indivíduo.\n")


## ----------------------------------------------------------------------------
## 2.3  Fragilidade gama semiparamétrica — survival::coxph com frailty()
## ----------------------------------------------------------------------------
##
## NOTA: {frailtypack} oferece mais opções (risco de base penalizado,
##   fragilidade correlacionada), mas requer instalação separada.
##   O pacote {survival} (já carregado) implementa fragilidade gama e
##   log-normal via frailty() diretamente no coxph(), o que é suficiente
##   para o laboratório.
##
## Para dados com cluster natural (famílias, hospitais):
##   frailty(id_cluster, distribution = "gamma")
##
## Agrupamos em 500 famílias de 10 indivíduos para simular dependência.

dados_cluster <- dados %>%
  mutate(familia = as.factor(rep(1:500, each = 10)))

cat("\n=== 2.3 Fragilidade gama — coxph() com cluster de família ===\n")
cox_frail_gama <- coxph(
  Surv(tempo, evento) ~ X + frailty(familia, distribution = "gamma"),
  data = dados_cluster
)
print(summary(cox_frail_gama))

## Extrair θ (variância da fragilidade gama)
hist_gama   <- cox_frail_gama$history[[1]]$history
theta_gama  <- hist_gama[nrow(hist_gama), "theta"]
cat("\nθ gama (variância da fragilidade):", round(theta_gama, 4), "\n")
cat("Valor verdadeiro σ² da simulação:  ", sigma2, "\n")
cat("\nInterpretação: θ =", round(theta_gama, 4), "\n")
if (theta_gama < 0.01) {
  cat("  θ ≈ 0 — resultado esperado nesta simulação.\n")
  cat("  As 500 'famílias' foram formadas agrupando indivíduos independentes\n")
  cat("  em sequência: não há dependência real intra-cluster, então o modelo\n")
  cat("  corretamente estima θ ≈ 0 (sem heterogeneidade entre clusters além de X).\n")
  cat("  Para observar θ > 0: simular dados com Z compartilhado por família, ou\n")
  cat("  usar dados reais (pbc, veteran na Parte 2 — lab2_dados_reais.R).\n")
} else {
  cat("  θ > 0 → heterogeneidade significativa entre famílias além de X.\n")
  cat("  Cada família compartilha um Z que multiplica seu risco de base.\n")
  cat("  θ → 0 reduziria o modelo ao Cox padrão sem fragilidade.\n")
}


## ----------------------------------------------------------------------------
## 2.4  Fragilidade log-normal — survival::coxph com frailty()
## ----------------------------------------------------------------------------
##
## TEORIA (slide 20):
##   Z = exp(W),   W ~ N(0, σ²)
##   Sem forma fechada para a Laplaciana → integração numérica (EM penalizado)
##   Cauda mais pesada que a gama: captura melhor Z >> 1 (heterogeneidade extrema)

cat("\n=== 2.4 Fragilidade log-normal — coxph() ===\n")
cox_frail_logn <- coxph(
  Surv(tempo, evento) ~ X + frailty(familia, distribution = "gaussian"),
  # "gaussian" no survival = fragilidade log-normal (W ~ N(0, θ))
  data = dados_cluster
)
print(summary(cox_frail_logn))

hist_logn  <- cox_frail_logn$history[[1]]$history
theta_logn <- hist_logn[nrow(hist_logn), "theta"]
cat("\nθ log-normal:", round(theta_logn, 4), "\n")

cat("\nNota técnica: distribution='gaussian' no coxph() corresponde a fragilidade\n")
cat("  log-normal porque o efeito aleatório entra na escala log do hazard:\n")
cat("  µ(x|W) = exp(W) · µ₀(x),  W ~ N(0, θ)  ↔  Z = exp(W) ~ Log-Normal\n")


## ----------------------------------------------------------------------------
## 2.5  Comparando distribuições: gama vs. log-normal
## ----------------------------------------------------------------------------
##
## Comparamos as densidades para θ estimados em 2.3 e 2.4.
## A log-normal tem cauda direita mais pesada — captura indivíduos com Z >> 1.
## O AIC/BIC do coxph orienta a escolha entre as duas especificações.

cat("\n=== 2.5 Comparação gama vs. log-normal ===\n")
cat("AIC — gama:      ", round(AIC(cox_frail_gama), 2), "\n")
cat("AIC — log-normal:", round(AIC(cox_frail_logn), 2), "\n")
cat("Menor AIC indica melhor ajuste. Diferença pequena → as duas são similares.\n")

z_seq <- seq(0.01, 5, length.out = 500)

dist_gama <- dgamma(z_seq, shape = 1/theta_gama, rate = 1/theta_gama)
# Log-normal: W ~ N(0, θ_logn), Z = exp(W)
# E[Z] = exp(θ_logn/2) ≠ 1 exatamente, mas próximo para θ pequeno
# Para manter E[Z]=1: meanlog = -θ_logn/2
dist_logn <- dlnorm(z_seq, meanlog = -theta_logn/2, sdlog = sqrt(theta_logn))

fig2 <- data.frame(
  z      = rep(z_seq, 2),
  dens   = c(dist_gama, dist_logn),
  modelo = rep(c("Gama", "Log-normal"), each = 500)
) |>
  ggplot(aes(x = z, y = dens, color = modelo, linetype = modelo)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "grey50") +
  scale_color_manual(values = c("#E64B35", "#4DBBD5")) +
  coord_cartesian(xlim = c(0, 4), ylim = c(0, 2)) +
  labs(
    title    = "Figura 2 — Distribuições de fragilidade estimadas",
    subtitle = "Gama vs. Log-normal: mesma média (1), caudas diferentes",
    x        = "Z (fragilidade)", y = "Densidade",
    color    = "Distribuição", linetype = "Distribuição",
    caption  = paste0(
      "theta_gama = ", round(theta_gama, 3),
      "  |  theta_logn = ", round(theta_logn, 3),
      "\nLinhas verticais em Z=1 (média populacional)."
    )
  )
print(fig2)

## PERGUNTA PARA DISCUSSÃO:
##   → A log-normal tem cauda mais pesada. Em que contextos isso importa?
##   → O AIC favorece qual distribuição nos dados simulados?
##   → O que acontece com o AIC se aumentarmos σ² na simulação para 1.5?


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
