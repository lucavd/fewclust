# fewclust 0.4.0: un pacchetto R classico, in stile tidyverse

Preparato il 2026-11-02 su richiesta del PI: "rendilo un pacchetto R classico, non un ibrido inutilizzabile; molto facile da usare, default robusti, intuitivo, messaggi chiari e utili; stile tidyverse, compatibile con quei pacchetti; non toccare il paper". Si parte dalla 0.3.1 (commit `626393b`). Il ramo del manoscritto non si tocca.

**Stato di questo file.** Le sezioni fino a "Test" sono scritte e committate *prima* del codice, insieme ai test e agli script in `dev/checks_0.4.0/`. Gli esiti si aggiungono alla fine.

## Cosa cambia

1. **Una sola funzione di stima, con i dati per primi** (stile tidyverse, va nelle pipe):
   `fit_fewclust(data, formula, id, time = NULL, method = c("auto", "mmm", "gee"), corstr = c("ar1", "exchangeable", "independence"), correction = c("fg", "md", "kc", "none"), conf_level = 0.95, control = fewclust_control())`.
   `id` e `time` si danno come nomi di colonna, nudi o fra virgolette (tidyselect), anche dentro funzioni dell'utente con `{{ }}`.
   - `method = "auto"` (default): il modello marginalizzato multilivello con il controllo di convergenza e la t su N - p gradi di libertà; se non dà una stima usabile, la GEE con correlazione AR(1) e correzione di Fay e Graubard, con un messaggio. È la regola della 0.3.1, invariata.
   - `method = "mmm"` o `"gee"` per sceglierne uno; per la GEE `corstr` e `correction`.
   - Il riferimento è sempre la t su N - p gradi di libertà (la scelta normale/t si toglie: la t era il default e l'unica raccomandata).
2. **Oggetto e metodi standard.** Classe `"fewclust"` con `print()`, `summary()`, `coef()`, `vcov()`, `confint()`, `nobs()`, `formula()`, `predict()` (scala logit o probabilità, anche su `newdata`), e i metodi broom `tidy()`, `glance()`, `augment()`, che restituiscono tibble.
3. **Messaggi.** Nessun numero dello studio nei messaggi di esecuzione; ogni messaggio dice cosa succede e cosa fare.
   - **Errori** (`cli::cli_abort()`, classe `"fewclust_error"`), solo quando non si può stimare: esito non binario o costante, esito che varia in un solo gruppo, separazione, più coefficienti che gruppi, una sola osservazione per gruppo, tempi ripetuti nel gruppo, termini a effetti casuali nella formula, nessun metodo con una stima usabile (il messaggio suggerisce `flip_test()`).
   - **Avvisi** (`cli::cli_warn()`, classe `"fewclust_warning"`), quando il risultato è meno affidabile: meno di 10 gruppi (prima era un rifiuto); GEE con 2 occasioni e meno di 15 gruppi; correzione diversa da Fay e Graubard; coefficienti senza inferenza (stima o errore standard ≥ 10); modello marginalizzato non convergente; buchi nelle serie o tempi non equispaziati.
   - **Messaggi** (`cli::cli_inform()`), informativi: ricaduta sulla GEE; `time` assente.
   - **Note nella stampa**, non condizioni: il controllo di convergenza è intervenuto; 2-3 occasioni e al massimo 25 gruppi (potenza bassa, `flip_test()`); con circa 10 gruppi il modello marginalizzato tende a essere conservativo.
4. **Si tolgono** il motore di decisione esposto e gli oggetti "rifiuto": `auto_cluster_fit()`, `cluster_strategies()`, `cluster_diagnose()`, `compare_strategies()`, i grafici, le strategie del modello misto (Wald e rapporto di verosimiglianza; lme4 esce dalle dipendenze). Le misure del disegno restano interne e alimentano errori, avvisi e note.
5. **`flip_test(data, formula, id, time = NULL, term = NULL, reps = 1999, enum_max = 12, seed = NULL)`**: stesso calcolo, dati per primi, avviso se non calcolabile, metodo `tidy()`.
6. **Invariati**: `gee_vcov()`, il motore del modello marginalizzato, il controllo di convergenza (default `check_coherence = TRUE`, decisione del PI), `simulate_clustered_binary()` (ora restituisce un tibble), `fc_trial`, le costanti dello studio (usate ora solo nella documentazione).
7. **Documentazione**: README e due vignette nuove ("Getting started", "What the simulation study found", con i numeri letti dalle costanti), NEWS con la corrispondenza fra vecchie e nuove funzioni, stile tidyverse del codice (styler) e lintr.

## Controlli (cancelli)

Script in `dev/checks_0.4.0/`; riferimento: la 0.3.1 installata da `626393b` in `~/MMM-package/.tmp040/lib031`; libreria di sviluppo `~/R/fewclust_dev_lib`; al massimo 4 core.

1. **Stessi numeri della 0.3.1** (`c1_identity.R`). 60 dataset simulati (10, 15, 25 e 65 gruppi; 2, 3, 4 e 8 occasioni; metà con una covariata; generati una volta con la 0.3.1 e riletti), per `method = "auto"`, `"mmm"` e la GEE con 3 correlazioni × 4 correzioni (14 stime per dataset): stima, errore standard, statistica, gradi di libertà e p-value di ogni coefficiente identici (`identical()`) alla strategia corrispondente della 0.3.1; dove la 0.3.1 rifiuta, la 0.4.0 dà un errore `fewclust_error`; `flip_test()` dà lo stesso p-value con lo stesso seme. Passa se tutto coincide.
2. **`gee_vcov()` e costanti invariati**: `gee_vcov()` identico bit per bit sugli stessi 60 dataset (3 correlazioni, `scale.fix` vero e falso, 4 tipi); `R/constants.R` con lo stesso md5 della 0.3.1.
3. **Nessun numero dello studio nei messaggi** (`tests/testthat/test-messages.R`): su una batteria di casi che attiva ogni errore, avviso, messaggio e nota, nessun testo contiene "design point", "Table", "scenario", "2,376", "mechanism", "of the study".
4. **Compatibilità tidyverse** (`tests/testthat/test-tidy.R`): tibble in ingresso, pipe `|>`, nomi di colonna nudi e fra virgolette, `{{ }}` in una funzione dell'utente, `tidy()`/`glance()`/`augment()` restituiscono tibble con i nomi di broom, uso con `dplyr::group_by() |> dplyr::group_modify()` e con `purrr::map()`.
5. **Test del pacchetto**: tutti passano; quelli numerici del motore (`test-mmm.R`, `test-mmm-check.R`, `test-corrections.R`) con le stesse attese della 0.3.1.
6. **`R CMD check --as-cran`**: 0 errori, 0 warning.
7. **Stile**: `lintr::lint_package()` con i linter di default (stile tidyverse) non trova nulla in `R/` (escluso `R/constants.R`, generato).

## Test

Riscritti sulla nuova interfaccia: `test-fit.R`, `test-methods.R`, `test-messages.R`, `test-tidy.R`, `test-flip.R` (stesso riferimento dello studio), `test-gee.R`, `test-mmm.R` (solo le chiamate alla funzione principale), `test-separation.R` e `test-diagnose.R` (misure interne). Tolti: `test-decision.R`, `test-compare.R`, `test-plot.R`, `test-glmm.R`. Invariati: `test-mmm-check.R`, `test-corrections.R`, `test-simulate.R` (salvo la classe tibble).

## Esiti

Eseguiti il 2026-11-02 con `bash dev/checks_0.4.0/run_checks.sh` sul server (R 4.6.1, al massimo 4 core); uscite in `dev/checks_0.4.0/out/`.

| Cancello | Esito | Dettaglio |
|---|---|---|
| 1. Stessi numeri della 0.3.1 | PASS | 795 stime confrontate, 0 differenze (`identical()` su termine, stima, errore standard, statistica, gradi di libertà, p-value); i 45 rifiuti della 0.3.1 sono tutti errori `fewclust_error` nella 0.4.0 (14 "esito che varia in un solo gruppo", 2 modello marginalizzato senza stima usabile, 29 GEE senza stima usabile); `flip_test()` 60/60 identico (p-value e motivo) |
| 2. `gee_vcov()` e costanti | PASS | 60/60 dataset identici (3 correlazioni × `scale.fix` × 4 tipi); md5 di `R/constants.R` identico alla 0.3.1 |
| 3. Nessun numero dello studio nei messaggi | PASS | `test-messages.R` |
| 4. Compatibilità tidyverse | PASS | `test-tidy.R` |
| 5. Test | PASS | 13 file, 494 attese, 0 fallite, 0 errori; 2 avvisi attesi (vedi sotto) |
| 6. `R CMD check --as-cran` | PASS | 0 errori, 0 warning, 1 NOTE (validazione HTML del manuale saltata: HTML Tidy non installato sul server, come nella 0.3.1) |
| 7. Stile | PASS | 0 lint in `R/` (la 0.3.1 ne aveva 315) |

### Scostamenti dal piano, dichiarati

1. **Tre errori nei test pre-registrati, corretti dopo il primo giro** (errori del test, non del pacchetto):
   - `helper-data.R::diag_of()` non dava `NULL` di default a `group` e non passava l'ambiente del chiamante a `fc_parse()`: tre test di `test-diagnose.R` fallivano con "argument is missing" e "object 'g' not found". Corretto con `group = NULL` ed `env = parent.frame()`.
   - In `test-flip.R` tre test scrivevano `ft <- expect_warning(flip_test(...))`; in testthat 3 `expect_warning()` restituisce la condizione, non il valore. Riscritti come `expect_warning(ft <- flip_test(...))`; attese invariate.
   - In `test-methods.R` (predict) `unname()` era applicato a un solo lato del confronto. Corretto applicandolo a entrambi.
2. **`run_checks.sh`, passo lintr**: la libreria di sviluppo non era nel percorso, quindi lintr non vedeva il namespace del pacchetto e segnalava come "no visible global function" le funzioni interne (82 falsi positivi). Aggiunta `DEV` a `.libPaths()`.
3. **Cancello 7 con soppressioni mirate**: `object_name_linter` è disattivato (`# nolint start/end`) in `R/mmm.R`, `R/corrections.R` e `R/separation.R`, che usano la notazione matriciale degli articoli (`X`, `Z`, `W`, `V`, ...; i nomi `X`, `Z`, `W` sono anche argomenti usati dai test numerici del motore), e sugli argomenti di broom `conf.int`, `conf.level`, `type.predict` di `tidy()`/`augment()`. Nessun altro linter è disattivato. Lo stile è stato applicato con `styler::style_file()` (tidyverse) su tutti i file tranne `R/constants.R`.
4. **Avvisi nei test**: due test (`test-messages.R`, 10 gruppi; `test-separation.R`, modello senza intercetta) producono l'avviso "The optimiser of the marginalized multilevel model did not converge (nlm code 5/3)": è l'avviso previsto dal piano (modello marginalizzato non convergente), non un errore.
5. **Ricaduta con `corstr`/`correction` dell'utente**: con `method = "auto"` la ricaduta usa la GEE con `corstr` e `correction` dati (di default AR(1) e Fay-Graubard, come nella 0.3.1). Il piano non lo specificava.

### Cosa non è cambiato

Motori (modello marginalizzato, GEE, correzioni, sign-flip), controllo di convergenza con `check_coherence = TRUE` di default, criterio di usabilità, costanti dello studio. Il ramo del manoscritto (`~/MMM`) non è stato toccato.
