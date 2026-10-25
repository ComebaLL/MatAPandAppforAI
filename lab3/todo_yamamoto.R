#install.packages(c("caustests", "lubridate", "zoo"))
library(caustests)
library(lubridate)
library(zoo)

# Параметры
CSV_PATH <- "C:/забгу/магистратура/Матаппарат и его приложения для ИИ/lab3/energy_dataset.csv"
COLS     <- c("generation.wind.onshore", "total.load.actual", "price.actual")
PMAX     <- 15     # потолок для выбора p по информационному критерию
DMAX     <- 1      # число добавочных лагов
IC       <- 1      # 1 - AIC, 2 - BIC
NBOOT    <- 499    # число бутстрэп-повторов для p-value


# Данные: суточные средние в уровнях
load_data <- function(path, cols) {
  df <- read.csv(path, check.names = TRUE)        # пробелы в именах -> точки
  df$time <- ymd_hms(df$time, tz = "UTC")
  x <- as.data.frame(lapply(df[, cols], na.approx, na.rm = FALSE, rule = 2))
  agg <- aggregate(x, by = list(day = as.Date(df$time)), FUN = mean, na.rm = TRUE)
  na.omit(agg[, cols])
}


# Готовая функция: caustests, test = 1 (Тода-Ямамото)
# Пакет принимает пару рядов и проверяет оба направления сразу,
# поэтому прогоняем все пары.
run_caustests <- function(data, cols) {
  out <- list()
  pairs <- combn(cols, 2, simplify = FALSE)
  for (pr in pairs) {
    cat("\n=== Пара:", pr[1], "<->", pr[2], "===\n")
    res <- caustests(data[, pr], test = 1, pmax = PMAX, ic = IC,
                     dmax = DMAX, nboot = NBOOT, verbose = FALSE)
    print(res)
    out[[paste(pr, collapse = " | ")]] <- res
  }
  invisible(out)
}


main <- function() {
  set.seed(42)                                    # воспроизводимость бутстрэпа
  data <- load_data(CSV_PATH, COLS)
  cat("Наблюдений:", nrow(data), "\n")

  run_caustests(data, COLS)
}

main()