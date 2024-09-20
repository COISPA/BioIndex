#' Estimation of ALK
#'
#' @param te MEDITS or MEDITS-like TE table
#' @param sp species rubin code (MEDITS format, e.g. "MERLMER")
#' @param GEAR type of gear reported in the corresponding TA file
#' @param GSA reference GSA for the analysis
#' @param country reference country
#' @param years number of years to be considered in the analysis
#' @export

ALKf <- function(te, sp, GEAR, GSA, country = NA, years = 5) {

    if (FALSE) {
        TE <- te #[te$AGE <30, ] # read.table("D:\\OneDrive - Coispa Tecnologia & Ricerca S.C.A.R.L\\BLACK SEA\\2023 - Data preparation - RAPA\\-DATA-\\BGR\\Spring 2022\\TE.csv", sep = ";", header = TRUE)
        # TE <- read.table("D:\\Documents and Settings\\Utente\\Documenti\\__ DATI MEDITS AGGIORNATI __\\_____MEDITS DATA\\GSA18\\TE 2012-20 con letture GSA18.csv", sep = ";", header = TRUE)
        # TE$SEX <- "N"
        sp <- "RAPAVEN"
        GSA <- 29
        country <- "TUE"
        years <- 5
        GEAR <- "BEAMT" # "TBB", "GOC73", "TRAWL"

        sp <- "MULLSUR"
        GSA <- 18
        country <- NA
        years <- 5
        GEAR <- "GOC73" # "BEAMT", "TBB", "GOC73", "TRAWL"


        ALK(TE, GEAR = "BEAMT", GSA = 29, country = "BGR", sp = "RAPAVEN", years = 15)

        ALK(TE, GEAR = "GOC73", GSA = 18, country = NA, sp = "MULLBAR", years = 15)

        ALK(TE, GEAR = "GOC73", GSA = 18, country = NA, sp = "MULLSUR", years = 5)

        ALK(TE, GEAR = "GOC73", GSA = 18, country = NA, sp = "MERLMER", years = 15)
    }

    TE <- te

    vbf3 <- function(t, Linf, k, t0) {
        Linf * (1 - exp(-k * (t - t0)))
    }

    wd <- getwd()
    sspp <- sp
    genus <- substr(unique(sspp), 1, 4)
    species <- substr(unique(sspp), 5, 7)

    TE$SEX <- toupper(TE$SEX)

    TE$COUNTRY <- as.character(TE$COUNTRY)
    if (all(is.na(country))) {
        country <- as.character(unique(TE$COUNTRY))
        TE <- TE[TE$AREA == GSA & TE$COUNTRY %in% country, ]
    } else {
        TE <- TE[TE$AREA == GSA & TE$COUNTRY %in% country, ]
    }
    data_C <- TE[TE$GENUS == genus & TE$SPECIES == species &
                     !is.na(TE$AGE) &
                     as.character(TE$AGE) != " " &
                     as.character(TE$AGE) != "" &
                     as.character(TE$AGE) != "UR" &
                     as.character(TE$AGE) != "NR" &
                     as.character(TE$AGE) != "-1" &
                     TE$SEX %in% c("F", "M", "N"), ]
    data_C$SEX <- "C"

    data <- TE[TE$GENUS == genus & TE$SPECIES == species &
                   !is.na(TE$AGE) &
                   as.character(TE$AGE) != " " &
                   as.character(TE$AGE) != "" &
                   as.character(TE$AGE) != "UR" &
                   as.character(TE$AGE) != "NR" &
                   as.character(TE$AGE) != "-1" &
                   TE$SEX %in% c("F", "M", "I", "N"), ]

    last_year <- max(TE$YEAR)
    years <- c((last_year - (years-1)):last_year)
    data <- rbind(data, data_C)
    data <- data[data$YEAR %in% years, ]
    data$AGE <- as.numeric(data$AGE)
    data <- data[!is.na(data$AGE), ]

    if (nrow(data) > 1) {
        sexes <- unique(data$SEX)
        # sexes <- sexes[sexes != "I"]
        s <- 3
        params <- list()
        for (s in 1:length(sexes)) {
            df <- data[data$SEX %in% c(sexes[s], "I"), ]
            if (nrow(df) > 0) {
                if (unique(GEAR)[1] %in% c("BEAMT", "TBB")) {
                    # ss <- suppressWarnings(vbStarts(SHELL_WIDTH~AGE,data=df,plot=FALSE))
                    ss <- list()
                    ss$Linf <- max(df$LENGTH_CLASS) / 0.95
                    if (length(df[df$AGE == 1, "LENGTH_CLASS"])>0){
                        ss$K <- mean(df[df$AGE == 1, "LENGTH_CLASS"]) / ss$Linf
                    } else {
                        ss$K <- 0.310 # (Sağlam, N., Sağlam, C., & Sağlam, Y. (2015))
                    }
                    ss$t0 <- -0.5
                    error <- try(mod3 <- nls(LENGTH_CLASS ~ vbf3(AGE, Linf, k, t0), start = c(Linf = abs(ss$Linf), k = abs(ss$K), t0 = ss$t0), data = df, algorithm = "default"),silent=TRUE)
                    if (!is(error, "try-error")) {
                        mod3 <- nls(LENGTH_CLASS ~ vbf3(AGE, Linf, k, t0), start = c(Linf = abs(ss$Linf), k = abs(ss$K), t0 = ss$t0), data = df, algorithm = "default")
                        sum <- summary(mod3)
                        sum
                        AIC(mod3)
                        Linf <- sum$parameters[1]
                        k <- sum$parameters[2]
                        t0 <- sum$parameters[3]
                        params[[s]] <- data.frame(species = sp, sex = sexes[s], Linf = as.numeric(round(Linf, 4)), k = as.numeric(round(k, 4)), t0 = as.numeric(round(t0, 4)), notes = "")
                        pred <- data.frame(AGE = seq(0, max(df$AGE), 0.1))
                        pred$LENGTH <- predict(mod3, newdata = data.frame(AGE = pred$AGE))
                        p <- ggplot() +
                            geom_point(data = df, aes(x = AGE, y = LENGTH_CLASS), colour = "red", stat = "identity") +
                            geom_line(data = pred, aes(x = AGE, y = LENGTH), colour = "blue",size=1) +
                            ggtitle(paste("ALK", sexes[s], sspp, sep = " - ")) +
                            labs(col = "Year") +
                            xlab("Age") +
                            theme(legend.position = "none") +
                            ylab(paste("Length (mm)"))
                        print(p)
                        ggsave(file.path(wd, "output","ALK", paste0("ALK_", sspp, "_", sexes[s], ".jpg")), dpi = 300, width = 9, height = 9, units = "in")
                    } else {
                        params[[s]] <- data.frame(species = sp, sex = sexes[s], Linf = NA, k = NA, t0 = NA, notes = "model not converged")
                        p <- ggplot() +
                            geom_point(data = df, aes(x = AGE, y = LENGTH_CLASS), colour = "red", stat = "identity") +
                            ggtitle(paste("ALK", sexes[s], sspp, sep = " - ")) +
                            labs(col = "Year") +
                            xlab("Age") +
                            theme(legend.position = "none") +
                            ylab(paste("Length (mm)"))
                        print(p)
                        ggsave(file.path(wd, "output","ALK", paste0("ALK_", sspp, "_", sexes[s], ".jpg")), dpi = 300, width = 9, height = 9, units = "in")
                    }
                }

                if (unique(GEAR)[1] %in% c("GOC73", "TRAWL")) {{ ss <- list()
                ss$Linf <- max(df$LENGTH_CLASS) / 0.95
                ss$K <- mean(df[df$AGE == 1, "LENGTH_CLASS"]) / ss$Linf
                ss$t0 <- -0.5
                # ss <- suppressWarnings(vbStarts(LENGTH_CLASS~AGE,data=df,plot=TRUE))

                error <- try(mod3 <- nls(LENGTH_CLASS ~ vbf3(AGE, Linf, k, t0), start = c(Linf = abs(ss$Linf), k = abs(ss$K), t0 = ss$t0), data = df, algorithm = "default"), silent = TRUE)
                if (!is(error, "try-error")) {
                    mod3 <- nls(LENGTH_CLASS ~ vbf3(AGE, Linf, k, t0), start = c(Linf = abs(ss$Linf), k = abs(ss$K), t0 = ss$t0), data = df, algorithm = "default")
                    sum <- summary(mod3)
                    sum
                    AIC(mod3)
                    Linf <- sum$parameters[1]
                    k <- sum$parameters[2]
                    t0 <- sum$parameters[3]
                    params[[s]] <- data.frame(species = sp, sex = sexes[s], Linf = as.numeric(round(Linf, 4)), k = as.numeric(round(k, 4)), t0 = as.numeric(round(t0, 4)), notes = "")
                    pred <- data.frame(AGE = seq(0, max(df$AGE), 0.1))
                    pred$LENGTH <- predict(mod3, newdata = data.frame(AGE = pred$AGE))
                    p <- ggplot() +
                        geom_point(data = df, aes(x = AGE, y = LENGTH_CLASS), colour = "red", stat = "identity") +
                        geom_line(data = pred, aes(x = AGE, y = LENGTH), colour = "blue",size=1) +
                        ggtitle(paste("ALK", sexes[s], sspp, sep = " - ")) +
                        labs(col = "Year") +
                        xlab("Age") +
                        theme(legend.position = "none") +
                        ylab(paste("Length (mm)"))
                    print(p)
                    ggsave(file.path(wd, "output","ALK", paste0("ALK_", sspp, "_", sexes[s], ".jpg")), dpi = 300, width = 9, height = 9, units = "in")
                } else {
                    params[[s]] <- data.frame(species = sp, sex = sexes[s], Linf = NA, k = NA, t0 = NA, notes = "model not converged")
                    p <- ggplot() +
                        geom_point(data = df, aes(x = AGE, y = LENGTH_CLASS), colour = "red", stat = "identity") +
                        ggtitle(paste("ALK", sexes[s], sspp, sep = " - ")) +
                        labs(col = "Year") +
                        xlab("Age") +
                        theme(legend.position = "none") +
                        ylab(paste("Length (mm)"))
                    print(p)
                    ggsave(file.path(wd, "output","ALK", paste0("ALK_", sspp, "_", sexes[s], ".jpg")), dpi = 300, width = 9, height = 9, units = "in")
                } }}
            } else {
                cat(paste0("Not age data for '", sexes[s], "' sex"))
            }
        }
        df.params <- do.call(rbind, params)
        write.table(df.params, file.path(wd, "output","ALK", paste0("ALK_", sspp, "_summary_table.csv")), sep=";", row.names=FALSE)
        return(df.params)
    } else {
        message("Not enougth data to plot ALK")
    }
} ### close function ALK
