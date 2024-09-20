#' Estimation of LW relationship
#'
#' @param TE MEDITS or MEDITS-like TE table
#' @param sp species rubin code (MEDITS format, e.g. "MERLMER")
#' @param GEAR type of gear reported in the corresponding TA file
#' @param GSA reference GSA for the analysis
#' @param country reference country
#' @param n_records minimum number of records to perform the analysis
#' @export
LWf <- function(TE, sp, GEAR, GSA, country=NA, n_records = 10) {
    if (FALSE) {
        sp <- "MERLMER"
        GEAR <- "GOC73"
        country <- country_analysis
        n_records <- 5



        TE=TE
        sp=sspp
        GEAR=GEAR
        GSA=GSA
        country = country_analysis
        years = n_years
        n_records=10

    }

    genus <- substr(unique(sspp), 1, 4)
    species <- substr(unique(sspp), 5, 7)

    sspp <- sp
    TE$COUNTRY <- as.character(TE$COUNTRY)
    if (all(is.na(country))) {
        country <- as.character(unique(TE$COUNTRY))
        TE <- TE[TE$AREA == GSA & TE$COUNTRY %in% country, ]
    } else {
        TE <- TE[TE$AREA == GSA & TE$COUNTRY %in% country, ]
    }

    if (unique(GEAR)[1] %in% c("BEAMT", "TBB")) {
        { #####  INDIVIDUAL BRUT WEIGHT - start #####
            data_iw <- TE[TE$GENUS == genus & TE$SPECIES == species &
                              as.character(TE$INDIVIDUAL_BRUT_WEIGHT) != "ND", ]
            data_iw$INDIVIDUAL_BRUT_WEIGHT <- as.numeric(data_iw$INDIVIDUAL_BRUT_WEIGHT)

            data <- data_iw
            colnames(data)[which(colnames(data) == "INDIVIDUAL_BRUT_WEIGHT")] <- "INDIVIDUAL_WEIGHT"
            var <- "Individual brut weight"

            if (nrow(data) > 1) {
                sexes <- unique(data$SEX)

                ###########   LW directory creation ###########################
                wd <- getwd()
                setwd(paste(wd, ("output/"), sep = "/"))
                dir_LW <- "LW"

                if (file.exists(dir_LW)) {
                    setwd(wd)
                } else {
                    dir.create(file.path(getwd(), dir_LW), showWarnings = FALSE)
                    setwd(wd)
                }
                LW_table_tot <- list()
                sex <- "M"
                for (sex in sexes) {
                    years <- unique(data$YEAR)
                    spdf <- data.frame(matrix(NA, ncol = 2, nrow = 1))
                    colnames(spdf) <- c("genus", "species")
                    spdf[, "genus"] <- substr(sspp, 1, 4)
                    spdf[, "species"] <- substr(sspp, 5, 7)

                    LW_table <- data.frame(matrix(nrow = length(years), ncol = 5))
                    colnames(LW_table) <- c("species", "sex", "a", "b", "year")
                    LW_table$species <- sspp
                    LW_table$sex <- sex


                    for (year in years) {


                        ###########   sex directory creation ###########################
                        setwd(paste(wd, ("output/LW/"), sep = "/"))
                        dir_sex <- sex

                        if (file.exists(dir_sex)) {
                            setwd(wd)
                        } else {
                            dir.create(file.path(getwd(), dir_sex), showWarnings = FALSE)
                            setwd(wd)
                        }
                        #################################################################
                        datay <- data[data$SEX == sex & data$YEAR == year & data$GENUS == genus & data$SPECIES == species, ]

                        # LC definition --------------------------------------------
                        LC_range <- c(min(datay$LENGTH_CLASS), max(datay$LENGTH_CLASS))

                        class_code <- unique(datay[datay$GENUS != -1, "LENGTH_CLASSES_CODE"])
                        if (class_code == "m") {
                            cl.size <- 1
                        } else {
                            if (class_code == "0") {
                                cl.size <- 5
                            } else {
                                if (class_code == "1") {
                                    cl.size <- 10
                                }
                            }
                        }

                        df <- data.frame(X = datay$LENGTH_CLASS, Y = datay$INDIVIDUAL_WEIGHT)
                        if (nrow(df) > n_records) {
                            model.lm <- lm(log(Y) ~ log(X), data = df)
                            a <- model.lm$coefficients[1]
                            b <- model.lm$coefficients[2]

                            if (class(try(nls(Y ~ a * (X^b), data = df, start = list(a = exp(a), b = b)),silent=TRUE)) != "try-error") {
                                mod <- nls(Y ~ a * (X^b), data = df, start = list(a = exp(a), b = b))
                                lc_list <- seq(LC_range[1], LC_range[2], cl.size)
                                prediction <- data.frame(matrix(nrow = length(lc_list), ncol = 2))
                                colnames(prediction) <- c("X", "pred")
                                prediction$X <- lc_list
                                prediction$pred <- predict(mod, newdata = prediction)

                                jpeg(file = paste(getwd(), "/output/LW/", sex, "/", sspp, "-", sex, "-", year, "-LW_Brut_weight.jpg", sep = ""), width = 20, height = 15, bg = "white", units = "cm", res = 200)
                                plot(df, pch = 16, cex = .5, main = paste(sspp, sex, sep = " - "), xlab = "length (mm)", ylab = "weight (g)")
                                lines(prediction$X, prediction$pred, col = "red", lwd = 2)
                                legend("topleft", c(paste("a = ", round(coef(mod)[[1]], 6)), paste("b = ", round(coef(mod)[[2]], 4))))
                                dev.off()

                                LW_table[which(years == year), 3] <- coef(mod)[[1]]
                                LW_table[which(years == year), 4] <- coef(mod)[[2]]
                                LW_table[which(years == year), 5] <- year
                            } else {
                                LW_table[which(years == year), 3] <- NA
                                LW_table[which(years == year), 4] <- NA
                                LW_table[which(years == year), 5] <- year
                                warning(paste("Not enougth data to estimate the length-weigth parameters for ", sspp, " - sex = ", sex, "\n", sep = ""))
                            } # try mod
                        } else {
                            LW_table[which(years == year), 3] <- NA
                            LW_table[which(years == year), 4] <- NA
                            LW_table[which(years == year), 5] <- year
                            warning(paste("Not enougth data to estimate the length-weigth parameters for ", sspp, " - sex = ", sex, "\n", sep = ""))
                        } # nrow di df > di n_records
                    } # years

                    LW_table_tot[[sex]] <- LW_table
                    write.table(LW_table, paste(getwd(), "/output/LW/", sex, "/LW_", sspp, "_", sex, "_", "Brut_Weight", ".csv", sep = ""), sep = ";", row.names = FALSE)
                }
                LW_table_final_brut <- do.call(rbind, LW_table_tot)
                LW_table_final_brut$notes <- "Brut weight"
                write.table(LW_table_final_brut, paste(getwd(), "/output/LW/LW_", sspp, "_", "Brut_Weight", ".csv", sep = ""), sep = ";", row.names = FALSE)
            } else {
                message("Not enougth data for length-weight analysis")
            }
        } #####  INDIVIDUAL BRUT WEIGHT - end #####


        { #####  CLEAN WEIGHT - start #####

            if (class(TE$CLEAN_WEIGHT)=="character"){
                data_cw <- TE[TE$GENUS == genus & TE$SPECIES == species &
                                  as.character(TE$CLEAN_WEIGHT != "ND"), ]
                data_cw$CLEAN_WEIGHT <- as.numeric(data_cw$CLEAN_WEIGHT)
            } else if(class(TE$CLEAN_WEIGHT)=="numeric") {
                data_cw <- TE[TE$GENUS == genus & TE$SPECIES == species &
                                  !is.na(TE$CLEAN_WEIGHT), ]
                data_cw$CLEAN_WEIGHT <- as.numeric(data_cw$CLEAN_WEIGHT)
            }


            data <- data_cw
            colnames(data)[which(colnames(data) == "CLEAN_WEIGHT")] <- "INDIVIDUAL_WEIGHT"
            var <- "Clean weight"

            if (nrow(data) > 1) {
                sexes <- unique(data$SEX)

                ###########   LW directory creation ###########################
                wd <- getwd()
                setwd(paste(wd, ("output/"), sep = "/"))
                dir_LW <- "LW"

                if (file.exists(dir_LW)) {
                    setwd(wd)
                } else {
                    dir.create(file.path(getwd(), dir_LW), showWarnings = FALSE)
                    setwd(wd)
                }
                sex="F"
                for (sex in sexes) {
                    years <- unique(data$YEAR)
                    spdf <- data.frame(matrix(NA, ncol = 2, nrow = 1))
                    colnames(spdf) <- c("genus", "species")
                    spdf[, "genus"] <- substr(sspp, 1, 4)
                    spdf[, "species"] <- substr(sspp, 5, 7)

                    LW_table <- data.frame(matrix(nrow = length(years), ncol = 5))
                    colnames(LW_table) <- c("species", "sex", "a", "b", "year")
                    LW_table$species <- sspp
                    LW_table$sex <- sex

                    year <- years[1]
                    for (year in years) {


                        ###########   sex directory creation ###########################
                        setwd(paste(wd, ("output/LW/"), sep = "/"))
                        dir_sex <- sex

                        if (file.exists(dir_sex)) {
                            setwd(wd)
                        } else {
                            dir.create(file.path(getwd(), dir_sex), showWarnings = FALSE)
                            setwd(wd)
                        }
                        #################################################################
                        datay <- data[data$SEX == sex & data$YEAR == year & data$GENUS == genus & data$SPECIES == species, ]

                        # LC definition --------------------------------------------
                        LC_range <- c(min(datay$LENGTH_CLASS), max(datay$LENGTH_CLASS))

                        class_code <- unique(datay[datay$GENUS != -1, "LENGTH_CLASSES_CODE"])
                        if (class_code == "m") {
                            cl.size <- 1
                        } else {
                            if (class_code == "0") {
                                cl.size <- 5
                            } else {
                                if (class_code == "1") {
                                    cl.size <- 10
                                }
                            }
                        }

                        df <- data.frame(X = datay$LENGTH_CLASS, Y = datay$INDIVIDUAL_WEIGHT)
                        if (nrow(df) > n_records) {
                            model.lm <- lm(log(Y) ~ log(X), data = df)
                            a <- model.lm$coefficients[1]
                            b <- model.lm$coefficients[2]

                            if (class(try(nls(Y ~ a * (X^b), data = df, start = list(a = exp(a), b = b)),silent=TRUE)) != "try-error") {
                                mod <- nls(Y ~ a * (X^b), data = df, start = list(a = exp(a), b = b))

                                lc_list <- seq(LC_range[1], LC_range[2], cl.size)
                                prediction <- data.frame(matrix(nrow = length(lc_list), ncol = 2))
                                colnames(prediction) <- c("X", "pred")
                                prediction$X <- lc_list
                                prediction$pred <- predict(mod, newdata = prediction)

                                jpeg(file = paste(getwd(), "/output/LW/", sex, "/", sspp, "-", sex, "-", year, "-LW_Clean_Weight.jpg", sep = ""), width = 20, height = 15, bg = "white", units = "cm", res = 200)
                                plot(df, pch = 16, cex = .5, main = paste(sspp, sex, sep = " - "), xlab = "length (mm)", ylab = "weight (g)")
                                lines(prediction$X, prediction$pred, col = "red", lwd = 2)
                                legend("topleft", c(paste("a = ", round(coef(mod)[[1]], 6)), paste("b = ", round(coef(mod)[[2]], 4))))
                                dev.off()

                                LW_table[which(years == year), 3] <- coef(mod)[[1]]
                                LW_table[which(years == year), 4] <- coef(mod)[[2]]
                                LW_table[which(years == year), 5] <- year
                            } else {
                                LW_table[which(years == year), 3] <- NA
                                LW_table[which(years == year), 4] <- NA
                                LW_table[which(years == year), 5] <- year
                                warning(paste("Not enougth data to estimate the length-weigth parameters for ", sspp, " - sex = ", sex, "\n", sep = ""))
                            } # try mod
                        } else {
                            LW_table[which(years == year), 3] <- NA
                            LW_table[which(years == year), 4] <- NA
                            LW_table[which(years == year), 5] <- year
                            warning(paste("Not enougth data to estimate the length-weigth parameters for ", sspp, " - sex = ", sex, "\n", sep = ""))
                        } # nrow di df > di n_records
                    } # years

                    LW_table_tot[[sex]] <- LW_table
                    write.table(LW_table, paste(getwd(), "/output/LW/", sex, "/LW_", sspp, "_", sex, "_", "Clean_Weight", ".csv", sep = ""), sep = ";", row.names = FALSE)
                }
                LW_table_final_clean <- do.call(rbind, LW_table_tot)
                LW_table_final_clean$notes <- "Clean weight"
                write.table(LW_table_final_clean, paste(getwd(), "/output/LW/LW_", sspp, "_", "Clean_Weight", ".csv", sep = ""), sep = ";", row.names = FALSE)

                if (exists("LW_table_final_brut") & exists("LW_table_final_clean")){
                    table_final <- rbind(LW_table_final_brut,LW_table_final_clean)
                } else if (!(exists("LW_table_final_brut") & exists("LW_table_final_clean"))) {
                    table_final <- LW_table_final_clean
                } else if ((exists("LW_table_final_brut") & !exists("LW_table_final_clean"))) {
                    table_final <- LW_table_final_brut
                } else {
                    table_final <- NULL
                }
                write.table(table_final, paste(getwd(), "/output/LW/LW_", sspp, ".csv", sep = ""), sep = ";", row.names = FALSE)
                return(LW_table)
            } else {
                message("Not enougth data for length-weight analysis")
            }

        } #### CLEAN WEIGHT - end
    } else if (unique(GEAR)[1] %in% c("GOC73", "TRAWL")) {
        data <- TE[TE$GENUS == genus & TE$SPECIES == species &
                       as.character(TE$INDIVIDUAL_WEIGHT) != "ND", ]

        if (nrow(data) > 1) {
            sexes <- unique(data$SEX)

            ###########   LW directory creation ###########################
            wd <- getwd()
            setwd(paste(wd, ("output/"), sep = "/"))
            dir_LW <- "LW"

            if (file.exists(dir_LW)) {
                setwd(wd)
            } else {
                dir.create(file.path(getwd(), dir_LW), showWarnings = FALSE)
                setwd(wd)
            }
            sex="I"
            for (sex in sexes) {
                years <- unique(data$YEAR)
                spdf <- data.frame(matrix(NA, ncol = 2, nrow = 1))
                colnames(spdf) <- c("genus", "species")
                spdf[, "genus"] <- substr(sspp, 1, 4)
                spdf[, "species"] <- substr(sspp, 5, 7)

                LW_table <- data.frame(matrix(nrow = length(years), ncol = 5))
                colnames(LW_table) <- c("species", "sex", "a", "b", "year")
                LW_table$species <- sspp
                LW_table$sex <- sex
                year <- 2012
                for (year in years) {

                    ###########   sex directory creation ###########################
                    setwd(paste(wd, ("output/LW/"), sep = "/"))
                    dir_sex <- sex

                    if (file.exists(dir_sex)) {
                        setwd(wd)
                    } else {
                        dir.create(file.path(getwd(), dir_sex), showWarnings = FALSE)
                        setwd(wd)
                    }
                    #################################################################
                    datay <- data[data$SEX == sex & data$YEAR == year & data$GENUS == genus & data$SPECIES == species, ]

                    # LC definition --------------------------------------------
                    LC_range <- c(min(datay$LENGTH_CLASS), max(datay$LENGTH_CLASS))

                    class_code <- unique(datay[datay$GENUS != -1, "LENGTH_CLASSES_CODE"])
                    if (class_code == "m") {
                        cl.size <- 1
                    } else {
                        if (class_code == "0") {
                            cl.size <- 5
                        } else {
                            if (class_code == "1") {
                                cl.size <- 10
                            }
                        }
                    }

                    df <- data.frame(X = datay$LENGTH_CLASS, Y = datay$INDIVIDUAL_WEIGHT)
                    if (nrow(df) > n_records) {
                        model.lm <- lm(log(Y) ~ log(X), data = df)
                        a <- model.lm$coefficients[1]
                        b <- model.lm$coefficients[2]

                        if (class(try(nls(Y ~ a * (X^b), data = df, start = list(a = exp(a), b = b)),silent=TRUE)) != "try-error") {
                            mod <- nls(Y ~ a * (X^b), data = df, start = list(a = exp(a), b = b))

                            lc_list <- seq(LC_range[1], LC_range[2], cl.size)
                            prediction <- data.frame(matrix(nrow = length(lc_list), ncol = 2))
                            colnames(prediction) <- c("X", "pred")
                            prediction$X <- lc_list
                            prediction$pred <- predict(mod, newdata = prediction)

                            jpeg(file = paste(getwd(), "/output/LW/", sex, "/", sspp, "-", sex, "-", year, "-LW.jpg", sep = ""), width = 20, height = 15, bg = "white", units = "cm", res = 200)
                            plot(df, pch = 16, cex = .5, main = paste(sspp, sex, sep = " - "), xlab = "length (mm)", ylab = "weight (g)")
                            lines(prediction$X, prediction$pred, col = "red", lwd = 2)
                            legend("topleft", c(paste("a = ", round(coef(mod)[[1]], 6)), paste("b = ", round(coef(mod)[[2]], 4))))
                            dev.off()

                            LW_table[which(years == year), 3] <- coef(mod)[[1]]
                            LW_table[which(years == year), 4] <- coef(mod)[[2]]
                            LW_table[which(years == year), 5] <- year
                        } else {
                            LW_table[which(years == year), 3] <- NA
                            LW_table[which(years == year), 4] <- NA
                            LW_table[which(years == year), 5] <- year
                            warning(paste("Not enougth data to estimate the length-weigth parameters for ", sspp, " - sex = ", sex, "\n", sep = ""))
                        } # try mod
                    } else {
                        LW_table[which(years == year), 3] <- NA
                        LW_table[which(years == year), 4] <- NA
                        LW_table[which(years == year), 5] <- year
                        warning(paste("Not enougth data to estimate the length-weigth parameters for ", sspp, " - sex = ", sex, "\n", sep = ""))
                    } # nrow di df > di n_records
                } # years
                write.table(LW_table, paste(getwd(), "/output/LW/", sex, "/LW_", sspp, "_", sex, "_", "Individual_Weight", ".csv", sep = ""), sep = ";", row.names = FALSE)
            } # sexes
        } else {
            message("Not enougth data for length-weight analysis")
        }
    }
} ### close function LW
