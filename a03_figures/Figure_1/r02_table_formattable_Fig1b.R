library(formattable)
library(readr)
table_1 <- read_delim("d02a_table_genome_info.txt", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)
table_1 <- table_1 %>% drop_na(ID)
table_1 <- table_1 %>% 
  mutate_if(is.character, ~replace_na(.,""))

unit.scale = function(x) (x - min(x)) / (max(x) - min(x))


formattable(table_1, 
            align = c("l","r", "r", "l", "l"),
            list(`ID` = formatter("span", style = ~ style(color = "grey", font.weight = "bold")),
                 `host` = formatter("span", style = ~ style(color = "black", font.style = "italic")),
                 `genome size (bp)` = color_bar("#C0E9C8", fun = unit.scale),
                 `# of contigs` = color_bar("#C0E9C8", fun = unit.scale),
                 `N50` = color_bar("#C0E9C8", fun = unit.scale), 
                 `GC content (%)` = color_tile("#caf4ff", "#00b7f1"), 
                 `Repetitive elements (%)` = color_tile("#caf4ff", "#00b7f1"), 
                 `gene count` = color_bar("#C0E9C8", fun = unit.scale), 
                 `BUSCO score (%)` = color_tile("#caf4ff", "#00b7f1"), 
                 `accessory genome (%)` = color_tile("#caf4ff", "#00b7f1"), 
                 `RIP (%)` = color_tile("#caf4ff", "#00b7f1")))

formattable(table_1, 
            align = c("l","r", "r", "l", "l", "l", "l", "l", "l", "l", "l", "l"),
            list(`ID` = formatter("span", style = ~ style(color = "grey", font.weight = "bold")),
                 `host` = formatter("span", style = ~ style(color = "black", font.style = "italic")),
                 `genome size (bp)` = color_bar("#C0E9C8", fun = unit.scale),
                 `# of contigs` = color_bar("#C0E9C8", fun = unit.scale),
                 `N50` = color_bar("#C0E9C8", fun = unit.scale), 
                 `GC content (%)` = color_bar("#87CEEB", fun = unit.scale), 
                 `Repetitive elements (%)` = color_bar("#87CEEB", fun = unit.scale), 
                 `gene count` = color_bar("#C0E9C8", fun = unit.scale), 
                 `BUSCO score (%)` = color_bar("#87CEEB", fun = unit.scale), 
                 `accessory genome (%)` = color_bar("#87CEEB", fun = unit.scale), 
                 `RIP (%)` = color_bar("#87CEEB", fun = unit.scale)))
