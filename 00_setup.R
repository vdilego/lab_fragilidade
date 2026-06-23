################################################################################
##  00_setup.R — Configuração do ambiente
##  Lab: Heterogeneidade, Fragilidade e Riscos Competitivos
##  Prof. Vanessa Di Lego · CEDEPLAR/UFMG · 2026
################################################################################

## ── Pacotes CRAN ─────────────────────────────────────────────────────────────
pkgs_cran <- c(
  # Análise de sobrevivência (fragilidade nativa, CIF, Fine & Gray)
  "survival",

  # Fragilidade semiparamétrica gama e log-normal (Rondeau et al.)
  "frailtypack",

  # Modelos paramétricos com fragilidade (Gompertz + gama)
  "flexsurv",

  # Fine & Gray clássico via crr() — alternativa ao finegray() do survival
  "cmprsk",

  # Manipulação e visualização
  "tidyverse",
  "ggplot2",
  "patchwork",   # combinar múltiplos gráficos
  "broom",       # tidy() para objetos survfit e coxph
  "ggsci",       # paletas de cores para publicação (jco, nejm, lancet)
  "scales"       # formatação de eixos (percent_format, comma)
)

# Instala apenas os pacotes ainda não instalados
pkgs_faltando <- pkgs_cran[!pkgs_cran %in% installed.packages()[, "Package"]]
if (length(pkgs_faltando) > 0) {
  install.packages(pkgs_faltando, dependencies = TRUE)
}

## ── Carregamento ──────────────────────────────────────────────────────────────
invisible(lapply(pkgs_cran, library, character.only = TRUE))

## ── Tema ggplot global ────────────────────────────────────────────────────────
theme_demog <- theme_bw(base_size = 13) +
  theme(
    legend.position   = "bottom",
    legend.title      = element_text(face = "bold"),
    strip.background  = element_rect(fill = "#2C3E50", color = NA),
    strip.text        = element_text(color = "white", face = "bold"),
    plot.title        = element_text(face = "bold", size = 13),
    plot.subtitle     = element_text(color = "grey40", size = 10),
    plot.caption      = element_text(color = "grey55", size = 8, hjust = 0),
    panel.grid.minor  = element_blank(),
    axis.text         = element_text(size = 10)
  )
theme_set(theme_demog)

## ── Paletas de cores padrão ───────────────────────────────────────────────────
# Consistentes com os slides do laboratório
COR_AZUL  <- "#2C3E50"
COR_TEAL  <- "#1A7F6E"
COR_CORAL <- "#C0392B"
COR_OURO  <- "#D4A843"
COR_CINZA <- "#6C757D"

## ── Reprodutibilidade ─────────────────────────────────────────────────────────
set.seed(2026)

cat("Ambiente configurado com sucesso!\n")
cat("Pacotes carregados:", paste(pkgs_cran, collapse=", "), "\n")
cat("Versão do R:", R.version$version.string, "\n")
cat("Versão do survival:", packageVersion("survival"), "\n")
