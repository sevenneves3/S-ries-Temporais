library(rmet)
library(tidyverse)
library(forecast)
library(tseries)

dados <- rmet::inmet_read(years = 2000:2026, stations = "A001")

glimpse(dados)

# Criando o data frame série mensal
serie_mensal <- dados %>%
  mutate(mes = lubridate::floor_date(datetime, unit = "month")) %>%
  group_by(mes) %>%
  summarise(temp_media = mean(temp_dry_c, na.rm = TRUE), .groups = "drop") %>%
  arrange(mes)

# Pegando o primeiro mês contabilizado
data_inicio <- min(serie_mensal$mes)

# Converte o vetor de médias mensais num objeto ts (série temporal)
ts_temp <- ts(serie_mensal$temp_media,
              start = c(lubridate::year(data_inicio), lubridate::month(data_inicio)),
              frequency = 12)

ts_temp


# Dividindo os dados para treino e teste
n_total <- length(ts_temp)
h_teste <- 12

ts_treino <- window(ts_temp, end = time(ts_temp)[n_total - h_teste])
ts_teste <- window(ts_temp, start = time(ts_temp)[n_total - h_teste + 1])


# Análise exploratória
autoplot(ts_temp) +
  labs(title = "Temperatura média mensal - Brasília (A001)",
       x = "Ano", y = "Temperatura média (°C)") +
  theme_minimal()

ggsubseriesplot(ts_treino) +
  labs(title = "Subséries mensais de temperatura", y = "Temperatura média (°C)") +
  theme_minimal()

serie_mensal %>%
  mutate(mes_nome = lubridate::month(mes, label = TRUE, abbr = TRUE)) %>%
  ggplot(aes(x = mes_nome, y = temp_media)) +
  geom_boxplot(fill = "steelblue", alpha = 0.6) +
  labs(title = "Distribuição mensal da temperatura média", x = "Mês", y = "Temperatura média (°C)") +
  theme_minimal()

medias_por_mes <- serie_mensal %>%
  mutate(mes_num = lubridate::month(mes)) %>%
  group_by(mes_num) %>%
  summarise(media = mean(temp_media, na.rm = TRUE)) %>%
  arrange(media)

medias_por_mes


#Decomposição sazonal
stl_decomp <- stl(ts_treino, s.window = "periodic")
autoplot(stl_decomp) +
  labs(title = "Decomposição STL da série de treino") +
  theme_minimal()

# ACF/PACF da série original
ggAcf(ts_treino, lag.max = 48) + labs(title = "ACF - série original") + theme_minimal()
ggPacf(ts_treino, lag.max = 48) + labs(title = "PACF - série original") + theme_minimal()


# Testes de estacionariedade
adf_result <- tseries::adf.test(ts_treino)
kpss_result <- tseries::kpss.test(ts_treino)

adf_result
kpss_result

nsdiffs(ts_treino)
ndiffs(ts_treino)

ts_treino_diff <- diff(ts_treino, lag = 12)

ggAcf(ts_treino_diff, lag.max = 48) + labs(title = "ACF - série com diferença sazonal") + theme_minimal()
ggPacf(ts_treino_diff, lag.max = 48) + labs(title = "PACF - série com diferença sazonal") + theme_minimal()


# Modelos candidatos
modelo_auto <- auto.arima(ts_treino, stepwise = FALSE, approximation = FALSE)

modelo_1 <- Arima(ts_treino, order = c(1, 0, 1), seasonal = list(order = c(0, 1, 2), period = 12), include.drift = TRUE)
modelo_2 <- Arima(ts_treino, order = c(1, 0, 0), seasonal = list(order = c(0, 1, 2), period = 12), include.drift = TRUE)
modelo_3 <- Arima(ts_treino, order = c(0, 0, 1), seasonal = list(order = c(0, 1, 1), period = 12), include.drift = TRUE)

comparacao <- tibble(
  modelo = c("SARIMA candidato 1", "SARIMA candidato 2", "SARIMA candidato 3", "auto.arima"),
  especificacao = c(
    paste0("(", paste(modelo_1$arma[c(1,6,2)], collapse = ","), ")(",
           paste(modelo_1$arma[c(3,7,4)], collapse = ","), ")[12]"),
    paste0("(", paste(modelo_2$arma[c(1,6,2)], collapse = ","), ")(",
           paste(modelo_2$arma[c(3,7,4)], collapse = ","), ")[12]"),
    paste0("(", paste(modelo_3$arma[c(1,6,2)], collapse = ","), ")(",
           paste(modelo_3$arma[c(3,7,4)], collapse = ","), ")[12]"),
    paste0("(", paste(modelo_auto$arma[c(1,6,2)], collapse = ","), ")(",
           paste(modelo_auto$arma[c(3,7,4)], collapse = ","), ")[12]")
  ),
  AIC = c(AIC(modelo_1), AIC(modelo_2), AIC(modelo_3), AIC(modelo_auto)),
  BIC = c(BIC(modelo_1), BIC(modelo_2), BIC(modelo_3), BIC(modelo_auto))
) %>%
  arrange(AIC)

comparacao

modelo_final <- modelo_1
summary(modelo_final)

# Diagnóstico de resíduos
checkresiduals(modelo_final)

ljung_box <- Box.test(residuals(modelo_final), lag = 24, type = "Ljung-Box", fitdf = length(modelo_final$coef))
ljung_box


# Validação no conjunto de teste
previsao_teste <- forecast(modelo_final, h = h_teste)

autoplot(previsao_teste) +
  autolayer(ts_teste, series = "Observado") +
  labs(title = "Previsão vs. valores observados - conjunto de teste",
       x = "Ano", y = "Temperatura média (°C)") +
  theme_minimal()

acuracia_teste <- accuracy(previsao_teste, ts_teste)
acuracia_teste

# Previsão final até dezembro de 2026
modelo_completo <- Arima(ts_temp, order = arimaorder(modelo_final)[1:3],
                          seasonal = list(order = arimaorder(modelo_final)[4:6], period = 12),
                          include.drift = TRUE)

ultimo_mes <- time(ts_temp)[length(ts_temp)]
horizonte_final <- round((2026 + 11/12 - ultimo_mes) * 12)

previsao_final <- forecast(modelo_completo, h = horizonte_final)

autoplot(previsao_final, include = 60) +
  labs(title = "Previsão da temperatura média mensal até dezembro de 2026",
       x = "Ano", y = "Temperatura média (°C)") +
  theme_minimal()

tabela_previsao <- tibble(
  mes = zoo::as.yearmon(time(previsao_final$mean)),
  previsao = as.numeric(previsao_final$mean),
  li_80 = previsao_final$lower[, 1],
  ls_80 = previsao_final$upper[, 1],
  li_95 = previsao_final$lower[, 2],
  ls_95 = previsao_final$upper[, 2]
)

tabela_previsao
