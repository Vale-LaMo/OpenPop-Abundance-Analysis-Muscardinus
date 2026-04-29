# supplementary mark tables

library(tidyverse)
library(knitr)
library(gt)

##---- Table S2 ----
# 1. Incolliamo i dati grezzi (puoi anche caricarli da un .csv se lo hai salvato)
raw_text <- "
 1:Phi:(Intercept)        2.9862211       0.5217797       1.9635328       4.0089093    
    2:Phi:Seasonhibernatio  -2.8681912       0.9402158      -4.7110142      -1.0253683    
    3:Phi:valleyS           -0.9364610       0.5742246      -2.0619413       0.1890192    
    4:Phi:Seasonhibernatio   1.3590916       1.6784464      -1.9306634       4.6488466    
    5:p:(Intercept)          0.3140275       1.4504659      -2.5288857       3.1569407    
    6:p:time1.492772667542  -0.8729022       1.2587901      -3.3401308       1.5943263    
    7:p:time2.018396846254  -1.9287997       1.4133515      -4.6989688       0.8413693    
    8:p:time2.511169513797  -1.8049940       1.4309208      -4.6095988       0.9996108    
    9:p:time3.003942181340  -1.7214035       1.4427329      -4.5491600       1.1063530    
   10:p:time3.496714848883  -2.0524545       1.4651233      -4.9240962       0.8191873    
   11:p:time4.022339027595  -1.3728094       1.4553029      -4.2252030       1.4795843    
   12:p:time4.515111695137  -2.3800782       1.4866567      -5.2939254       0.5337690    
   13:p:time5.040735873850  -2.1519419       1.4821472      -5.0569505       0.7530667    
   14:p:time5.533508541392  -2.6268690       1.4982971      -5.5635315       0.3097934    
   15:p:time6.354796320630  -4.1560452       1.6215458      -7.3342750      -0.9778155    
   16:p:time13.31931668856  -1.3140609       1.4969399      -4.2480632       1.6199414    
   17:p:time14.33771353482  -1.8588956       1.5220421      -4.8420982       1.1243070    
   18:p:time15.32325886990  -2.1254832       1.5253786      -5.1152253       0.8642590    
   19:p:time16.01314060446  -1.6962543       1.5124844      -4.6607237       1.2682152    
   20:p:time16.50591327201  -1.9392277       1.5131703      -4.9050416       1.0265862    
   21:p:time17.03153745072  -1.2901901       1.5088432      -4.2475229       1.6671426    
   22:p:time17.52431011826  -3.3770700       1.5638024      -6.4421229      -0.3120172    
   23:p:time26.36136662286  -1.4330867       1.5493191      -4.4697521       1.6035788    
   24:p:time27.34691195795  -1.7953950       1.5681442      -4.8689577       1.2781677    
   25:p:time28.36530880420  -1.2280362       1.5509045      -4.2678091       1.8117367    
   26:p:time29.38370565045  -2.5230881       1.5584764      -5.5777020       0.5315258    
   27:p:time30.36925098554  -3.3831111       1.6072672      -6.5333548      -0.2328673    
   28:pent:(Intercept)      -0.1453132       0.6805308      -1.4791535       1.1885271    
   29:N:(Intercept)          2.2003419       0.4833715       1.2529338       3.1477501    
   30:N:groupS1             -0.6265363       0.8466219      -2.2859152       1.0328426    
   31:N:groupR2              0.1622083       0.6222202      -1.0573433       1.3817599    
   32:N:groupS2             -0.8773326       0.9334638      -2.7069217       0.9522565    
   33:N:groupR3              1.8651603       0.4910507       0.9027009       2.8276198    
   34:N:groupS3              0.9796527       0.5709910      -0.1394897       2.0987950"

# 2. Pulizia e formattazione
mark_table <- read.table(text = raw_text, col.names = c("ID_Raw", "Beta", "SE", "LCI", "UCI")) %>%
  separate(ID_Raw, into = c("ID", "Parameter", "Variable"), sep = ":", fill = "right") %>%
  mutate(
    # Arrotondamento a 3 decimali
    across(where(is.numeric), ~ round(., 3)),
    # Pulizia nomi variabili per renderli leggibili
    Variable = str_replace(Variable, "Seasonhibernatio", "Season (Winter)"),
    Variable = str_replace(Variable, "valleyS", "Valley (Valsavarenche)"),
    Variable = ifelse(is.na(Variable) | Variable == "(Intercept)", "Intercept", Variable)
  ) %>%
  select(Parameter, Variable, Beta, SE, LCI, UCI)

# 3. Visualizzazione (formato Markdown per il paper)
kable(mark_table, format = "markdown")

mark_table %>%
  gt() %>%
  # tab_header(
  #   title = md("**Model Comparison via LOO-CV**"),
  #   subtitle = "Sensitivity analysis of prior distributions"
  # ) %>%
  # Riduzione della dimensione del font e spaziatura
  tab_options(
    table.font.size = px(7),
    data_row.padding = px(2.5),
    column_labels.font.size = px(8),
    column_labels.font.weight = "bold"
  ) %>%
  fmt_number(
    columns = where(is.numeric),
    decimals = 3
  ) %>%
  cols_label(
    Parameter = "Parameter",
    Variable = "Variable",
    Beta = "Coefficient",
    SE = "SE",
    LCI = "LCI",
    UCI = "UCI"
  ) %>%
  # # Evidenzia il modello vincitore
  # tab_style(
  #   style = cell_fill(color = "lightgrey", alpha = 0.5),
  #   locations = cells_body(rows = 1)
  # ) %>%
  gtsave("outputs/POPAN/Table_S2.png")

## ---- Table S3 ----
raw_text <- "
  1:Phi:(Intercept)        2.5540876       0.5816921       1.4139710       3.6942042    
  2:Phi:Seasonhibernatio  -2.9426110       1.0028402      -4.9081779      -0.9770441    
  3:Phi:valleyS           -1.0546570       0.6192289      -2.2683457       0.1590317    
  4:Phi:altitude2          0.5467407       0.5660808      -0.5627776       1.6562591    
  5:Phi:altitude3          0.6616575       0.4556199      -0.2313576       1.5546726    
  6:Phi:Seasonhibernatio   1.3890082       1.6591799      -1.8629845       4.6410009    
  7:p:(Intercept)          0.3634801       1.4950129      -2.5667452       3.2937054    
  8:p:time1.492772667542  -0.9125517       1.3040929      -3.4685739       1.6434704    
  9:p:time2.018396846254  -1.9765477       1.4578965      -4.8340249       0.8809294    
 10:p:time2.511169513797  -1.8567709       1.4755430      -4.7488353       1.0352934    
 11:p:time3.003942181340  -1.7764322       1.4872222      -4.6913878       1.1385233    
 12:p:time3.496714848883  -2.1062945       1.5088833      -5.0637059       0.8511169    
 13:p:time4.022339027595  -1.4277326       1.4993311      -4.3664215       1.5109564    
 14:p:time4.515111695137  -2.4333170       1.5298494      -5.4318218       0.5651879    
 15:p:time5.040735873850  -2.2054869       1.5255081      -5.1954827       0.7845090    
 16:p:time5.533508541392  -2.6811360       1.5412676      -5.7020206       0.3397486    
 17:p:time6.354796320630  -4.2118990       1.6613486      -7.4681423      -0.9556557    
 18:p:time13.31931668856  -1.3985394       1.5392756      -4.4155197       1.6184409    
 19:p:time14.33771353482  -1.9469624       1.5658066      -5.0159433       1.1220186    
 20:p:time15.32325886990  -2.2133191       1.5689560      -5.2884730       0.8618347    
 21:p:time16.01314060446  -1.7860889       1.5563484      -4.8365319       1.2643542    
 22:p:time16.50591327201  -2.0262898       1.5568821      -5.0777787       1.0251991    
 23:p:time17.03153745072  -1.3844575       1.5523364      -4.4270368       1.6581219    
 24:p:time17.52431011826  -3.4604014       1.6058545      -6.6078763      -0.3129265    
 25:p:time26.36136662286  -1.5766442       1.5877873      -4.6887073       1.5354189    
 26:p:time27.34691195795  -1.9474373       1.6101288      -5.1032898       1.2084151    
 27:p:time28.36530880420  -1.3735129       1.5932863      -4.4963540       1.7493282    
 28:p:time29.38370565045  -2.6416807       1.6017028      -5.7810183       0.4976569    
 29:p:time30.36925098554  -3.4953366       1.6489992      -6.7273752      -0.2632980    
 30:pent:(Intercept)      -0.1172529       0.6839151      -1.4577265       1.2232208    
 31:N:(Intercept)          2.4531132       0.4787207       1.5148207       3.3914057    
 32:N:groupS1             -0.5651281       0.7795838      -2.0931125       0.9628562    
 33:N:groupR2             -0.0876417       0.6502300      -1.3620925       1.1868090    
 34:N:groupS2             -1.0834933       0.9420538      -2.9299188       0.7629322    
 35:N:groupR3              1.5644833       0.5000432       0.5843985       2.5445681    
 36:N:groupS3              0.7111074       0.5711303      -0.4083080       1.8305228"

 mark_table <- read.table(text = raw_text, col.names = c("ID_Raw", "Beta", "SE", "LCI", "UCI")) %>%
  separate(ID_Raw, into = c("ID", "Parameter", "Variable"), sep = ":", fill = "right") %>%
  mutate(
    # Arrotondamento a 3 decimali
    across(where(is.numeric), ~ round(., 3)),
    # Pulizia nomi variabili per renderli leggibili
    Variable = str_replace(Variable, "Seasonhibernatio", "Season (Winter)"),
    Variable = str_replace(Variable, "valleyS", "Valley (Valsavarenche)"),
    Variable = ifelse(is.na(Variable) | Variable == "(Intercept)", "Intercept", Variable)
  ) %>%
  select(Parameter, Variable, Beta, SE, LCI, UCI)

mark_table %>%
  gt() %>%
  tab_options(
    table.font.size = px(7),
    data_row.padding = px(2.5),
    column_labels.font.size = px(8),
    column_labels.font.weight = "bold"
  ) %>%
  fmt_number(
    columns = where(is.numeric),
    decimals = 3
  ) %>%
  cols_label(
    Parameter = "Parameter",
    Variable = "Variable",
    Beta = "Coefficient",
    SE = "SE",
    LCI = "LCI",
    UCI = "UCI"
  ) %>%
  gtsave("outputs/POPAN/Table_S3.png")

## ---- Table S4 ----
raw_text <- "
1:Phi:(Intercept)        1.9272503       0.4694153       1.0071962       2.8473044    
    2:Phi:Seasonhibernatio  -1.0150211       2.1897426      -5.3069167       3.2768745    
    3:Phi:altitude2          2.1242860       2.1129557      -2.0171071       6.2656792    
    4:Phi:altitude3          0.7292939       0.5711691      -0.3901976       1.8487854    
    5:Phi:Seasonhibernatio  -4.1367579       3.2038673      -10.416338       2.1428221    
    6:Phi:Seasonhibernatio  -1.1296498       2.2496462      -5.5389563       3.2796568    
    7:p:(Intercept)          0.1932945       1.3485526      -2.4498686       2.8364577    
    8:p:time1.492772667542  -0.8002734       1.1613064      -3.0764340       1.4758872    
    9:p:time2.018396846254  -1.8357879       1.3136196      -4.4104824       0.7389066    
   10:p:time2.511169513797  -1.7067888       1.3296689      -4.3129400       0.8993624    
   11:p:time3.003942181340  -1.6221988       1.3410416      -4.2506404       1.0062428    
   12:p:time3.496714848883  -1.9452709       1.3640133      -4.6187371       0.7281952    
   13:p:time4.022339027595  -1.2623131       1.3526648      -3.9135362       1.3889100    
   14:p:time4.515111695137  -2.2630762       1.3861360      -4.9799029       0.4537505    
   15:p:time5.040735873850  -2.0317620       1.3810458      -4.7386118       0.6750877    
   16:p:time5.533508541392  -2.5044455       1.3983208      -5.2451543       0.2362633    
   17:p:time6.354796320630  -4.0335020       1.5296879      -7.0316903      -1.0353136    
   18:p:time13.31931668856  -1.1902846       1.3949794      -3.9244443       1.5438752    
   19:p:time14.33771353482  -1.7402015       1.4242108      -4.5316548       1.0512518    
   20:p:time15.32325886990  -1.9942018       1.4280302      -4.7931412       0.8047375    
   21:p:time16.01314060446  -1.5566280       1.4142964      -4.3286491       1.2153931    
   22:p:time16.50591327201  -1.7967715       1.4149405      -4.5700549       0.9765118    
   23:p:time17.03153745072  -1.1432959       1.4103312      -3.9075451       1.6209534    
   24:p:time17.52431011826  -3.2353636       1.4688205      -6.1142519      -0.3564754    
   25:p:time26.36136662286  -1.2305931       1.4534070      -4.0792710       1.6180848    
   26:p:time27.34691195795  -1.6630371       1.4722524      -4.5486518       1.2225777    
   27:p:time28.36530880420  -1.0854599       1.4538697      -3.9350445       1.7641247    
   28:p:time29.38370565045  -2.3821718       1.4620931      -5.2478742       0.4835307    
   29:p:time30.36925098554  -3.2390123       1.5142217      -6.2068868      -0.2711378    
   30:pent:(Intercept)      -0.2051227       0.6717049      -1.5216644       1.1114189    
   31:N:(Intercept)          2.5241493       0.4722533       1.5985328       3.4497659    
   32:N:groupS1             -0.9351550       0.8025475      -2.5081481       0.6378381    
   33:N:groupR2             -0.1834332       0.6511308      -1.4596496       1.0927831    
   34:N:groupS2             -1.5711049       1.0415072      -3.6124591       0.4702493    
   35:N:groupR3              1.5820786       0.4942960       0.6132584       2.5508988    
   36:N:groupS3              0.3949237       0.5617593      -0.7061245       1.4959719"

 mark_table <- read.table(text = raw_text, col.names = c("ID_Raw", "Beta", "SE", "LCI", "UCI")) %>%
  separate(ID_Raw, into = c("ID", "Parameter", "Variable"), sep = ":", fill = "right") %>%
  mutate(
    # Arrotondamento a 3 decimali
    across(where(is.numeric), ~ round(., 3)),
    # Pulizia nomi variabili per renderli leggibili
    Variable = str_replace(Variable, "Seasonhibernatio", "Season (Winter)"),
    Variable = str_replace(Variable, "valleyS", "Valley (Valsavarenche)"),
    Variable = ifelse(is.na(Variable) | Variable == "(Intercept)", "Intercept", Variable)
  ) %>%
  select(Parameter, Variable, Beta, SE, LCI, UCI)

mark_table %>%
  gt() %>%
  tab_options(
    table.font.size = px(7),
    data_row.padding = px(2.5),
    column_labels.font.size = px(8),
    column_labels.font.weight = "bold"
  ) %>%
  fmt_number(
    columns = where(is.numeric),
    decimals = 3
  ) %>%
  cols_label(
    Parameter = "Parameter",
    Variable = "Variable",
    Beta = "Coefficient",
    SE = "SE",
    LCI = "LCI",
    UCI = "UCI"
  ) %>%
  gtsave("outputs/POPAN/Table_S4.png")
