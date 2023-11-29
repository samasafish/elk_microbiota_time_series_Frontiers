# functions for phyloseq
require(phyloseq);
require(tidyverse);
require(viridis);


#### plotting in phyloseq ####
color_abundance_plot = function (prevalence_df, physeq, tax_level= "Order", legend= FALSE){
  prevalence_df2 = subset(prevalence_df, prevalence_df[[tax_level]] %in% get_taxa_unique(physeq, tax_level));
  if(legend==FALSE){
    ggplot(prevalence_df2, aes(TotalAbundance, Prevalence / nsamples(physeq),color=prevalence_df2[[tax_level]])) +
      geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  
      geom_point(size = 2, alpha = 0.7) +
      scale_x_log10(labels= scales::comma) +  
      xlab("Total Abundance") + 
      ylab("Prevalence [Frac. Samples]") +
      facet_wrap(~c(tax_level)) + 
      theme(legend.position="none") #margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  }
  else{
    ggplot(prevalence_df2, aes(TotalAbundance, Prevalence / nsamples(physeq),color=prevalence_df2[[tax_level]])) +
      geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  
      geom_point(size = 2, alpha = 0.7) +
      scale_x_log10(labels= scales::comma) +  
      ylab("Prevalence [Frac. Samples]") +
      facet_wrap(~c(tax_level)) + 
      theme(legend.text = element_text(size = 6, margin = margin(r = 0, unit = 'pt')), 
            legend.key=element_blank(),legend.title = element_blank(),
            legend.position="bottom",legend.spacing.x = unit(-.1, 'cm'),
            legend.margin = margin(-5,0,0,-40)) #margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  }
}
#if error resize window or save as object
#add pallete option
color_abundance_plot2 = function (prevalence_df, physeq, tax_level= "", pallet_builder=TRUE, legend= FALSE, title= ""){
  #set color pallete
  require(ggplot2);require(phyloseq);require(wesanderson);
  if(pallet_builder==T){
    #if is not NULL
    pal_COLORS <- as.vector(wesanderson::wes_palette(length(phyloseq::get_taxa_unique(physeq, tax_level)), 
                                                     name = "GrandBudapest1", type = "continuous"))
    
    names(pal_COLORS) <- levels(prevalence_df[[tax_level]])
    
    low_letter <- tolower(substr(tax_level, start = 1, stop = 1))
    pal_COLORS[which(names(pal_COLORS)==paste0(low_letter, "__Unknown"))] <- "green" 
    pal_COLORS[which(names(pal_COLORS)=="NA")] <- "green" 
    colScale <- scale_colour_manual(name = tax_level, values = pal_COLORS)
  }
  else{colScale = NULL }
  
  if(legend==FALSE){
    ggplot(prevalence_df, aes(TotalAbundance, Prevalence / nsamples(physeq),color=prevalence_df[[tax_level]])) +
      colScale + ggtitle(title) +
      geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  
      geom_point(size = 2, alpha = 0.7) +
      scale_x_log10(labels= scales::comma) +  ylim(0,NA) +
      xlab("Total Abundance") + 
      ylab("Prevalence [Frac. Samples]") +
      facet_wrap(~c(tax_level)) + 
      theme(legend.position="none") #margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  }
  else{
    ggplot(prevalence_df, aes(TotalAbundance, Prevalence / nsamples(physeq), color=prevalence_df[[tax_level]])) +
      colScale + ggtitle(title) +
      geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  
      geom_point(size = 2, alpha = 0.7) +
      scale_x_log10(labels= scales::comma) +   ylim(0,NA) +
      ylab("Prevalence [Frac. Samples]") +
      facet_wrap(~c(tax_level)) + 
      theme(legend.text = element_text(size = 6, margin = margin(r = 0, unit = 'pt')), 
            legend.key=element_blank(),legend.title = element_blank(),
            legend.position="bottom",legend.spacing.x = unit(-.1, 'cm'),
            legend.margin = margin(-5,0,0,-40)) #margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
  }
}

# make data-frame of nreads for both ASVs and Samples and plot
plot_nreads_summary = function(physeq, log10y= F){ 
  require(phyloseq);require(ggplot2);
  title = "Total number of reads (log10)"
  readsums.df = data.frame(nreads = sort(taxa_sums(physeq), TRUE), sorted = 1:ntaxa(physeq), type = "ASVs")
  readsums.df = rbind(readsums.df, data.frame(nreads = sort(sample_sums(physeq),TRUE), sorted = 1:nsamples(physeq), type = "Samples"))
  if(log10y){
    ggplot(readsums.df, aes(x = sorted, y = nreads)) +
      geom_bar(fill = "Black", stat = "identity", alpha = 1) +
      geom_bar(data=readsums.df, stat = "identity",fill = "red", alpha = .5) +
      ggtitle(title) + facet_wrap(~type, 1, scales = "free") +
      scale_y_log10(labels= scales::comma)
  }
  else{
    title = "Total number of reads"
    ggplot(readsums.df, aes(x = sorted, y = nreads)) +
      geom_bar(fill = "Black", stat = "identity", alpha = 1) +
      geom_bar(data=readsums.df, stat = "identity",fill = "red", alpha = .5) +
      ggtitle(title) + facet_wrap(~type, 1, scales = "free") +
      scale_y_continuous(trans ="identity", labels= scales::comma)
  }
}

#compare and plot read depths for ASVs and Samples post transformation or filter
plot_compare_nreads_summary = function(physeq, physeq2, log10y= F, title=""){ 
  require(phyloseq);require(ggplot2);
  title2 = " (log10)"
  #original
  readsums.df = data.frame(nreads = sort(taxa_sums(physeq), TRUE), sorted = 1:ntaxa(physeq), type = "ASVs")
  readsums.df = rbind(readsums.df, data.frame(nreads = sort(sample_sums(physeq),TRUE), sorted = 1:nsamples(physeq), type = "Samples"))
  #filtered
  f.readsums.df = data.frame(nreads = sort(taxa_sums(physeq2), TRUE), sorted = 1:ntaxa(physeq2), type = "ASVs")
  f.readsums.df = rbind(f.readsums.df, data.frame(nreads = sort(sample_sums(physeq2),TRUE), sorted = 1:nsamples(physeq2), type = "Samples"))
  #plot
  if(log10y){
    ggplot(readsums.df, aes(x = sorted, y = nreads)) +
      geom_bar(fill = "Black", stat = "identity", alpha = 1) +
      geom_bar(data=f.readsums.df, stat = "identity",fill = "red", alpha = .5) +
      ggtitle(paste0(title, title2)) + 
      facet_wrap(~type, 1, scales = "free") +
      scale_y_log10(labels= scales::comma)
  }
  else{
    title = "Total number of reads after Variance Transformation"
    ggplot(readsums.df, aes(x = sorted, y = nreads)) + 
      geom_bar(fill = "Black", stat = "identity", alpha = 1) + 
      geom_bar(data=f.readsums.df,stat = "identity",fill = "red", alpha = .5) +
      ggtitle(paste0(title, title2)) + 
      facet_wrap(~type, 1, scales = "free") + 
      scale_x_continuous(labels = scales::comma)
  }
}

# make df to simplify counting specific taxonomic levels
prevalence_df = function(physeq){
  prevdf = apply(X = otu_table(physeq),
                 MARGIN = ifelse(taxa_are_rows(physeq), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(physeq),
                      tax_table(physeq))
  return(prevdf)
} 

#pies from phyloseq object
plot_pie_phyloseq = function (physeq, title = "", color_cat= "Species", legend = T){
  require(viridis);require(phyloseq);
  if(legend){
    plot_bar(physeq,x ="Grazer_Species", fill = color_cat) + 
      geom_bar(width=1, stat="identity") +
      coord_polar("y") +
      scale_fill_viridis(discrete=TRUE) +
      ggtitle(title) +
      theme_void() +
      theme(plot.title = element_text(hjust=0.5),
            legend.text = element_text(size = 7, margin = margin(r = .5, unit = 'pt')), 
            legend.key=element_blank(),legend.title = element_blank(),
            legend.position="right",legend.spacing.x = unit(0, 'cm'))
    # legend.margin = margin(0,0,0,0)) 
    #margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  } 
  else {
    plot_bar(physeq,x ="Grazer_Species", fill = color_cat) + 
      geom_bar(width=1, stat="identity") +
      coord_polar("y") +
      scale_fill_viridis(discrete=TRUE) +
      ggtitle(title) +
      theme_void() +
      theme(plot.title = element_text(hjust=0.5), legend.position="none")
  }
}

#pies from data.frame object (after ps.melt())
plot_pie_df = function(melted.df, title= ""){ 
  require(viridis);require(ggplot2);
  ggplot(data.frame(melted.df), 
         aes(x="Species", y=Abundance, fill=Species)) + 
    geom_bar(width=1, stat="identity") +
    coord_polar("y") +
    scale_fill_viridis(discrete=TRUE) +
    ggtitle(title) +
    theme_void() +
    theme(plot.title = element_text(hjust=0.5), 
          legend.title=element_blank())
}

#for collapsing low abundance into a single group (may need adapting depending on column names)
collapse_low_abundance = function(physeq, percent_collapse = 1){
  require(phyloseq);
  #   if(glom_rank=T) {
  #      glom.ps <- tax_glom(physeq, taxrank = glom_rank)
  #      physeq <- glom.ps
  #      data_pie<- psmelt(physeq); # create dataframe from phyloseq object
  #      data_pie$glom_rank <- as.character(data_pie$glom_rank); #convert to character
  #      data_pie$glom_rank[data_pie$glom_rank < percent_collapse] <- sprintf("< %d perc.abund.", percent_collapse);
  #      temp_tax_vec <- list(unique(data_pie$glom_rank));
  #      data_pie$glom_rank <- factor(data_pie$glom_rank, levels =unlist(temp_tax_vec));
  # }else if {
  data_pie = psmelt(physeq) # create dataframe from phyloseq object
  data_pie$Species <- as.character(data_pie$Species) #convert to character
  #simple way to rename phyla with < 1% abundance
  data_pie$Species[data_pie$Abundance < percent_collapse] <- sprintf("< %d perc. abund.", percent_collapse)
  temp_spp_vec <- list(unique(data_pie$Species))
  data_pie$Species <- factor(data_pie$Species, levels =unlist(temp_spp_vec))
  return(data_pie)
} #species only
collapse_low_abundance_levels = function(physeq, tax_level="Species", taxa_glom = T, percent_collapse = 1){
  require(phyloseq);
  if (taxa_glom==T) {
    #if the data is transformed into proportions this will work.. if not the percent collapse is actually counts
    physeq <- phyloseq::tax_glom(physeq, taxrank = tax_level)
    data_pie = phyloseq::psmelt(physeq) # create dataframe from phyloseq object
    data_pie[[tax_level]] <- as.character(data_pie[[tax_level]]) #convert to character
    #simple way to rename phyla with < 1% abundance
    data_pie[[tax_level]][data_pie$Abundance < percent_collapse] <- sprintf("< %d perc. abund.", percent_collapse)
    temp_spp_vec <- list(unique(data_pie[[tax_level]]))
    data_pie[[tax_level]] <- factor(data_pie[[tax_level]], levels =unlist(temp_spp_vec))
    return(data_pie)
  }
  else{
    data_pie = psmelt(physeq) # create dataframe from phyloseq object
    data_pie[[tax_level]] <- as.character(data_pie[[tax_level]]) #convert to character
    #simple way to rename phyla with < 1% abundance
    data_pie[[tax_level]][data_pie$Abundance < percent_collapse] <- sprintf("< %d perc. abund.", percent_collapse)
    temp_spp_vec <- list(unique(data_pie[[tax_level]]))
    data_pie[[tax_level]] <- factor(data_pie[[tax_level]], levels =unlist(temp_spp_vec))
    return(data_pie)
  } 
}
# any taxonomic level
#### plotting helper functions ####
pallet_builder_physeq = function(physeq, tax_level= tax_level, unknownColor = "black"){
  pal_COLORS <- as.vector(wesanderson::wes_palette(length(phyloseq::get_taxa_unique(physeq, tax_level)), 
                                                   name = "GrandBudapest1", type = "continuous"))
  names(pal_COLORS) <-  get_taxa_unique(physeq, taxonomic.rank = tax_level)
  low_letter <- tolower(substr(tax_level, start = 1, stop = 1))
  pal_COLORS[which(names(pal_COLORS)==paste0(low_letter, "__Unknown"))] <- unknownColor 
  pal_COLORS[which(names(pal_COLORS)=="NA")] <- "green" 
  colScale <- scale_colour_manual(name = tax_level, values = pal_COLORS)
  return(pal_COLORS)
} # for phyloseq object
pallet_builder_melt = function(melt.df, tax_level= tax_level, unknownColor = "black"){
  
  pal_COLORS<- as.vector(wesanderson::wes_palette(length(unique(melt.df[[tax_level]])),
                                                  name = "GrandBudapest1", type = "continuous"))
  names(pal_COLORS_test) <-  unique(melt.df[[tax_level]])
  low_letter <- tolower(substr(tax_level, start = 1, stop = 1))
  pal_COLORS[which(names(pal_COLORS)==paste0(low_letter, "__Unknown"))] <- unknownColor 
  pal_COLORS[which(names(pal_COLORS)=="NA")] <- "green" 
  colScale <- scale_colour_manual(name = tax_level, values = pal_COLORS)
  return(pal_COLORS)
} # for dataframe (melted phyloseq)
reorder_col_levels = function(melt.df, tax_level1 = tax_level, orderByCol = "Abundance"){
  order_levels = unique(melt.df[[tax_level1]][order(-melt.df[[orderByCol]])])
  melt.df[[tax_level1]] <- factor(melt.df[[tax_level1]], levels =order_levels, ordered = T)
  return(melt.df)
}
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}
# big bar plot functions ####
# this plot is obsolete since it displays counts without regard for sample number (samples are summed within grazer groups)
plot_bar_with_colorSet = function(physeq, tax_level="", legend=F, reorder=F, title = "title",
                                  collapseAbundPerc = NULL, collapseColor = "green",
                                  taxa_glom = F, legendLoc= "right", Legendcols= 2, 
                                  rmUnkTax = F, unknownColor = "black"){
  require(phyloseq);require(ggplot2);require(wesanderson);
  # currently only supporting wesAnderson colors
  if (rmUnkTax==T){
    cat("rmUnkTax is on... removing unkown taxa from specified tax_level with glom function\n")
    low_letter <- tolower(substr(tax_level, start = 1, stop = 1))
    unknown_tax_name <- paste0(low_letter, "__Unknown")
    unknown_tax_name <- c(NA, "", " ", "\t",unknown_tax_name)
    physeq<- tax_glom(physeq, taxrank=tax_level, NArm=TRUE, bad_empty=unknown_tax_name)
  }
  else {physeq = physeq}
  ##### collapse low abundance taxa?
  if (!is.null(collapseAbundPerc)){ #if not NULL
    physeq.melt <- collapse_low_abundance_levels(physeq, tax_level = tax_level, percent_collapse = collapseAbundPerc, taxa_glom = taxa_glom)
    #this function is from functions script and needs to be sourced first
    #also data should be transformed to proportions or the collapseAbundPerc is actually counts
  } 
  else {
    if(taxa_glom == T){
      physeq <- phyloseq::tax_glom(physeq = physeq,taxrank = tax_level)
      physeq.melt <- psmelt(physeq)
    }
    else{
      physeq.melt <- psmelt(physeq)}
  }
  #get pallet everytime  
  pal_COLORS<- as.vector(wesanderson::wes_palette(length(unique(physeq.melt[[tax_level]])), 
                                                  name = "GrandBudapest1", type = "continuous"))
  names(pal_COLORS) <-  unique(physeq.melt[[tax_level]])
  low_letter <- tolower(substr(tax_level, start = 1, stop = 1))
  pal_COLORS[which(names(pal_COLORS)==paste0(low_letter, "__Unknown"))] <- unknownColor 
  pal_COLORS[which(names(pal_COLORS)=="NA")] <- "green" 
  pal_COLORS[which(names(pal_COLORS)== sprintf("< %d perc. abund.", collapseAbundPerc))] <- collapseColor 
  colScale <- scale_colour_manual(name = tax_level, values = pal_COLORS)
  ### reorder?
  if (reorder == T){
    # order_levels = unique(physeq.melt[[tax_level]][order(-physeq.melt$Abundance)])
    # physeq.melt[[tax_level]] <- factor(physeq.melt[[tax_level]], levels =order_levels, ordered = T)
    physeq.melt <- reorder_col_levels(melt.df = physeq.melt, tax_level = tax_level,orderByCol = "Abundance")
  } 
  else {order_level = NULL}
  #legend
  if (legend == T){
    ggplot(physeq.melt, 
           aes(x=Grazer_Species, y= Abundance, factor(physeq.melt[[tax_level]]), 
               fill=factor(physeq.melt[[tax_level]]))) + 
      ggtitle(title) +
      geom_bar(aes(),stat = "identity")+
      scale_fill_manual(values = pal_COLORS)+ guides(fill=guide_legend(ncol=Legendcols)) +
      theme(legend.text = element_text(size = 6, margin = margin(r = 0, unit = 'pt')), 
            legend.key=element_blank(),legend.title = element_blank(),
            legend.position=legendLoc,legend.spacing.x = unit(0, 'cm'), #distance from legend block
            legend.margin = margin(0,0,0,-10)) #distance from bars
  } 
  else {
    ggplot(physeq.melt, 
           aes(x=Grazer_Species, y= Abundance, factor(physeq.melt[[tax_level]]), 
               fill=factor(physeq.melt[[tax_level]]))) + 
      ggtitle(title) +
      geom_bar(aes(),stat = "identity")+
      scale_fill_manual(values = pal_COLORS)+ theme(legend.position="none")
    
  }
}
#less options
group_100bar_plot = function(physeq, group = "Grazer_Species", tax_level = "Species", nTaxa = NULL,
                             title= "title", collapseAbundPerc = NULL, collapseColor = "green", taxa_glom = F,
                             legend = T, Legendcols = 2, legendLoc = "right", facetGrid = NULL,
                             rmUnkTax = F, rmTaxaList = NULL, unknownColor = "black"){
  require(phyloseq);require(ggplot2);require(wesanderson);
  cat("Takes a phyloseq object and produces bar plots... \n")
  cat("with proportional bars for each sample group\n")
  # filter to nTaxa
  if (!is.null(rmTaxaList)){ #if not NULL
    cat("rmTaxaList is on... removing taxa\n")
    physeq <- pop_taxa(physeq, badTaxa = rmTaxaList)
  }
  else {physeq = physeq}
  
  if (rmUnkTax==T){
    cat("rmUnkTax is on... removing unkown taxa from specified tax_level with glom function\n")
    low_letter <- tolower(substr(tax_level, start = 1, stop = 1))
    unknown_tax_name <- paste0(low_letter, "__Unknown")
    unknown_tax_name <- c(NA, "", " ", "\t",unknown_tax_name)
    physeq<- tax_glom(physeq, taxrank=tax_level, NArm=TRUE, bad_empty=unknown_tax_name)
  }
  else {physeq = physeq}
  
  if (!is.null(nTaxa)){ #if not NULL
    cat("nTaxa is on... displaying top ", nTaxa, " taxa\n")
    ntaxaNames <- names(sort(taxa_sums(physeq), decreasing=TRUE))[1:nTaxa] #names
    physeq <- prune_taxa(ntaxaNames, physeq) # prune on names
  }
  else {physeq = physeq}
  # if filtering with nTaxa the community proportions may be misrepresented 
  #get transformation
  merge.ps = merge_samples(physeq, group)
  sample_data(merge.ps)[[group]] <- levels(sample_data(physeq)[[group]])
  # make a copy of table for proportions by group
  otu.table <- otu_table(merge.ps)
  # apply proportions by species
  # merge.ps <- apply(merge.ps, MARGIN = 2, function(x){100*x/sum(x)})
  # apply proportions by GROUP
  otu.table <- apply(otu.table, MARGIN = 1, function(x){100*x/sum(x)})
  otu_table(merge.ps) <- otu_table(otu.table, taxa_are_rows=T)
  
  if (!is.null(collapseAbundPerc)){ #if not NULL
    cat("collapseAbundPerc is on...\ncollapsing taxa with less than", collapseAbundPerc, "percent abundance into single group\n")
    merge.melt <- collapse_low_abundance_levels(merge.ps, tax_level = tax_level, percent_collapse = collapseAbundPerc, taxa_glom = taxa_glom)
    #this function is from functions script and needs to be sources first
  } 
  else {
    if(taxa_glom == T){
      merge.ps <- phyloseq::tax_glom(physeq = merge.ps,taxrank = tax_level, NArm = F)
      merge.melt <- psmelt(merge.ps)
    }
    else{
      merge.melt <- psmelt(merge.ps)}
  }
  
  #get pallet everytime  
  pal_COLORS<- as.vector(wesanderson::wes_palette(length(unique(merge.melt[[tax_level]])), 
                                                  name = "GrandBudapest1", type = "continuous"))
  names(pal_COLORS) <-  unique(merge.melt[[tax_level]])
  low_letter <- tolower(substr(tax_level, start = 1, stop = 1))
  pal_COLORS[which(names(pal_COLORS)==paste0(low_letter, "__Unknown"))] <- unknownColor 
  pal_COLORS[which(names(pal_COLORS)=="NA")] <- "green" 
  pal_COLORS[which(names(pal_COLORS)== sprintf("< %d perc. abund.", collapseAbundPerc))] <- collapseColor 
  colScale <- scale_colour_manual(name = tax_level, values = pal_COLORS)
  #draw plot
  if (legend == T){
    cat("plotting with legend\n")
    if (!is.null(facetGrid)){
      cat("plotting with facet_grid ~", facetGrid,"\n")
      
      ggplot(merge.melt, 
             aes(x=merge.melt[[group]], y= Abundance,
                 fill=factor(merge.melt[[tax_level]]))) +
        ggtitle(title) + geom_bar(aes(),stat = "identity")+ xlab(group) + 
        facet_grid(cols = vars(merge.melt[["Grazer_Species"]]), scales = "free", labeller = label_bquote(`facetGrid`)) +
        scale_fill_manual(values = pal_COLORS)+ guides(fill=guide_legend(ncol=Legendcols)) +
        theme(legend.text = element_text(size = 6, margin = margin(r = 0, unit = 'pt')), 
              legend.key=element_blank(),legend.title = element_blank(),
              legend.position=legendLoc,legend.spacing.x = unit(0, 'cm'), #distance from legend block
              legend.margin = margin(0,0,0,-10)) #distance from bars
    }
    else {
      ggplot(merge.melt, 
             aes(x=merge.melt[[group]], y= Abundance, 
                 fill=factor(merge.melt[[tax_level]]))) + 
        ggtitle(title) + geom_bar(aes(),stat = "identity")+ xlab(group) +
        scale_fill_manual(values = pal_COLORS)+ guides(fill=guide_legend(ncol=Legendcols)) +
        theme(legend.text = element_text(size = 6, margin = margin(r = 0, unit = 'pt')), 
              legend.key=element_blank(),legend.title = element_blank(),
              legend.position=legendLoc,legend.spacing.x = unit(0, 'cm'), #distance from legend block
              legend.margin = margin(0,0,0,-10)) #distance from bars
    }
  } 
  else {
    cat("plotting without legend\n")
    if (!is.null(facetGrid)){
      cat("plotting with facet_grid ~", facetGrid,"\n")
      ggplot(merge.melt, 
             aes(x=merge.melt[[group]], y= Abundance, 
                 fill=factor(merge.melt[[tax_level]]))) + 
        facet_grid(cols = vars(merge.melt[["Grazer_Species"]]), scales = "free", labeller = label_bquote(`facetGrid`)) +
        ggtitle(title) + xlab(group) +
        geom_bar(aes(),stat = "identity")+
        scale_fill_manual(values = pal_COLORS)+ theme(legend.position="none")
    }
    else {
      ggplot(merge.melt, 
             aes(x=merge.melt[[group]], y= Abundance, 
                 fill=factor(merge.melt[[tax_level]]))) + 
        ggtitle(title) + xlab(group) +
        geom_bar(aes(),stat = "identity")+
        scale_fill_manual(values = pal_COLORS)+ theme(legend.position="none")
    }
  }
}
#more options
#also works for SampleID bars, Grazer bars
group_100bar_plot2 = function(physeq, group = "Grazer_Species", tax_level = "Species", nTaxa = NULL,
                              title= "title", legend = T, Legendcols = 2, legendLoc = "right", 
                              collapseAbundPerc = NULL, collapseColor = "green", taxa_glom = F,
                              facetGrid = NULL, renameUnkTax =F, reorder =F,
                              rmUnkTax = F, rmTaxaList = NULL, unknownColor = "black", 
                              colorPal = "wes", barColSeperator=F, reColorTax = NULL, reColorList = NULL){
  require(phyloseq);require(ggplot2);require(wesanderson);
  cat("Produces proportional bars for each sample or group in a physeq object.\n")
  physeq_original <- physeq
  levels_group <- levels(sample_data(physeq_original)[[group]])
  # cat("this is levels_group", levels_group, "\n")
  # filter to nTaxa
  
  # after glom...
  # sample_data(T20.vst.merge.grazers)$Grazer_Species <- levels(sample_data(top20.vst.ps)$Grazer_Species)
  if (!is.null(rmTaxaList)){ #if not NULL
    cat("rmTaxaList is on... removing taxa\n")
    physeq <- pop_taxa(physeq, badTaxa = rmTaxaList)
  }
  else {physeq = physeq}
  if (rmUnkTax==T){
    cat("rmUnkTax is on... removing unkown taxa from specified tax_level with glom function\n")
    low_letter <- tolower(substr(tax_level, start = 1, stop = 1))
    unknown_tax_name <- paste0(low_letter, "__Unknown")
    unknown_tax_name <- c(NA, "", " ", "\t",unknown_tax_name)
    physeq<- tax_glom(physeq, taxrank=tax_level, NArm=TRUE, bad_empty=unknown_tax_name)
  }
  else {physeq = physeq}
  
  if (!is.null(nTaxa)){ #if not NULL
    cat("nTaxa is on... displaying top ", nTaxa, " taxa summed across all groups\n")
    ntaxaNames <- names(sort(taxa_sums(physeq), decreasing=TRUE))[1:nTaxa] #names
    physeq <- prune_taxa(ntaxaNames, physeq) # prune on names
  }
  else {physeq = physeq}
  # if filtering with nTaxa the community proportions may be misrepresented 
  #get transformation
  merge.ps = merge_samples(physeq, group)
  sample_data(merge.ps)[[group]] <- levels(sample_data(physeq)[[group]])
  # sample_data(merge.ps)[[tax_level]] <- levels(sample_data(merge.ps)[[tax_level]])
  if (!is.null(facetGrid)){
    # sample_data(merge.ps)[[facetGrid]] <- levels(sample_data(merge.ps)[[facetGrid]])
    sample_data(merge.ps)[[facetGrid]] <- factor(x = sample_data(physeq)[[facetGrid]], 
                                                 levels = levels(sample_data(physeq)[[facetGrid]]))
  }
  
  # make a copy of table for proportions by group
  otu.table <- otu_table(merge.ps)
  # apply proportions by species
  # merge.ps <- apply(merge.ps, MARGIN = 2, function(x){100*x/sum(x)})
  # apply proportions by GROUP
  otu.table <- apply(otu.table, MARGIN = 1, function(x){100*x/sum(x)})
  otu_table(merge.ps) <- otu_table(otu.table, taxa_are_rows=T)
  
  if (!is.null(collapseAbundPerc)){ #if not NULL
    cat("collapseAbundPerc is on...\ncollapsing taxa with less than", collapseAbundPerc, "percent abundance into single group\n")
    merge.melt <- collapse_low_abundance_levels(merge.ps, tax_level = tax_level, 
                                                percent_collapse = collapseAbundPerc, taxa_glom = taxa_glom)
    
    if (!is.null(facetGrid)){ #reset factor levels for facet group
      merge_levels <- levels(merge.melt[[facetGrid]])
      cat("levels of facet variable ", merge_levels, "\n")
    }
  } 
  else {
    if(taxa_glom == T){
      merge.ps <- phyloseq::tax_glom(physeq = merge.ps,taxrank = tax_level, NArm = F)
      merge.melt <- psmelt(merge.ps)
    }
    else{
      merge.melt <- psmelt(merge.ps)}
  }
  
  if (renameUnkTax==T){
    cat("renameUnkTax is on... converting unkowns specified tax_level to Family+OTU ID\n")
    #get tax_level into variable
    low_letter <- tolower(substr(tax_level, start = 1, stop = 1))
    unknown_tax_name <- paste0(low_letter, "__Unknown")
    #index the unknowns at the specific taxonomic level and the coresponding IDs
    unknown_index <- which(merge.melt[[tax_level]]== paste0(low_letter, "__Unknown"))
    
    #set options
    family_names = F
    class_names = T
    Raw_OTU_names = F
    
    if(family_names == TRUE){
      # add OTU name at each unknown entry to the factor levels
      modified_OTU_names <- ""
      modified_OTU_names <- gsub(pattern = "\\D+?(16s)\\D+" ,replacement =  "OTU_" , x = as.character(merge.melt$OTU[unknown_index]), perl = T)
      modified_OTU_names <- paste(merge.melt$Family[unknown_index], modified_OTU_names, sep="_")
      modified_OTU_names <- gsub(pattern = "f__Unknown" ,replacement =  "Family_unk" , x = modified_OTU_names, perl = T)
    } 
    if(class_names == TRUE){
      # add OTU name at each unknown entry to the factor levels
      modified_OTU_names <- ""
      modified_OTU_names <- gsub(pattern = "\\D+?(16s)\\D+" ,replacement =  "OTU_" , x = as.character(merge.melt$OTU[unknown_index]), perl = T)
      modified_OTU_names <- paste(merge.melt$Family[unknown_index], modified_OTU_names, sep="_")
      modified_OTU_names <- paste(merge.melt$Class[unknown_index], modified_OTU_names, sep="_")
      modified_OTU_names <- gsub(pattern = "f__Unknown" ,replacement =  "Family_unk" , x = modified_OTU_names, perl = T)
      modified_OTU_names <- gsub(pattern = "c__Unknown" ,replacement =  "Class_unk" , x = modified_OTU_names, perl = T)
      # print(modified_OTU_names)
    } 
    if(Raw_OTU_names == TRUE){
      modified_OTU_names <- merge.melt$OTU[which(merge.melt[[tax_level]]== paste0(low_letter, "__Unknown"))]
      # modified_OTU_names <- unique(modified_OTU_names)
      # add OTU name at each unknown entry to the factor levels
      # levels(merge.melt[[tax_level]]) <- c(levels(merge.melt[[tax_level]]), modified_OTU_names)
      #   # add OTU name at each unknown entry
      # merge.melt[[tax_level]][unknown_index] <- as.character(merge.melt$OTU[unknown_index])
      # known_col_length <- length(unique(merge.melt[[tax_level]]))
      # Unk_col_length <- length(modified_OTU_names)
    }
    #use the new OTU names to modify the data frame and levels
    levels(merge.melt[[tax_level]]) <- c(levels(merge.melt[[tax_level]]), unique(modified_OTU_names))
    # add OTU name at each unknown entry
    merge.melt[[tax_level]][unknown_index] <- modified_OTU_names
    # print(modified_OTU_names)
    #for col palette
    known_col_length <- length(unique(merge.melt[[tax_level]]))
    Unk_col_length <- length(modified_OTU_names)
  }
  # when rename is off
  else{ 
    cat("renameUnkTax is off... Unknown taxa remain unnamed\n")
    low_letter <- tolower(substr(tax_level, start = 1, stop = 1))
    unknown_tax_name <- paste0(low_letter, "__Unknown")
    #index the unknowns at the specific taxonomic level and the coresponding IDs
    unknown_index <- which(merge.melt[[tax_level]]== paste0(low_letter, "__Unknown"))
    modified_OTU_names <- merge.melt$OTU[which(merge.melt[[tax_level]]== paste0(low_letter, "__Unknown"))]
    modified_OTU_names <- unique(modified_OTU_names)
    #need a length thats based on how many unknowns are included
    #this will be fed to the color palette
    known_col_length <- length(unique(merge.melt[[tax_level]]))
    Unk_col_length <- length(modified_OTU_names)
  }
  #order by abundance
  if (reorder ==T){
    cat("Reorder taxa by total abundance (global)\n")
    # reorder_col_levels = function(melt.df, tax_level = tax_level, orderByCol = "Abundance")
    merge.melt <- reorder_col_levels(merge.melt)
  }
  
  if (colorPal == "wes"){
    pal_COLORS<- as.vector(wesanderson::wes_palette(known_col_length, 
                                                    name = "GrandBudapest1", type = "continuous"))
    names(pal_COLORS) <-  unique(merge.melt[[tax_level]])
    low_letter <- tolower(substr(tax_level, start = 1, stop = 1))
    pal_COLORS[which(names(pal_COLORS)==paste0(low_letter, "__Unknown"))] <- unknownColor 
    pal_COLORS[which(names(pal_COLORS)=="NA")] <- "green" 
    pal_COLORS[which(names(pal_COLORS)== sprintf("< %d perc. abund.", collapseAbundPerc))] <- collapseColor 
    colScale <- scale_colour_manual(name = tax_level, values = pal_COLORS)
  }
  
  if (colorPal == "rainbow"){
    pal_COLORS<- rainbow(n = known_col_length, v = 0.7)
    names(pal_COLORS) <-  unique(merge.melt[[tax_level]])
    low_letter <- tolower(substr(tax_level, start = 1, stop = 1))
    pal_COLORS[which(names(pal_COLORS)==paste0(low_letter, "__Unknown"))] <- unknownColor 
    pal_COLORS[which(names(pal_COLORS)=="NA")] <- "green" 
    pal_COLORS[which(names(pal_COLORS)== sprintf("< %d perc. abund.", collapseAbundPerc))] <- collapseColor 
    colScale <- scale_colour_manual(name = tax_level, values = pal_COLORS)
  }
  if (colorPal == "qual"){
    require(RColorBrewer);
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    qual_col_pals <- qual_col_pals[-c(4,5,7,8),]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    pal_COLORS <- col_vector[1:known_col_length]
    names(pal_COLORS) <-  unique(merge.melt[[tax_level]])
    low_letter <- tolower(substr(tax_level, start = 1, stop = 1))
    pal_COLORS[which(names(pal_COLORS)==paste0(low_letter, "__Unknown"))] <- unknownColor 
    pal_COLORS[which(names(pal_COLORS)=="NA")] <- "green" 
    pal_COLORS[which(names(pal_COLORS)== sprintf("< %d perc. abund.", collapseAbundPerc))] <- collapseColor 
    colScale <- scale_colour_manual(name = tax_level, values = pal_COLORS)
  }
  
  if (colorPal == "wes2"){
    pal_COLORS<- as.vector(wesanderson::wes_palette(known_col_length, 
                                                    name = "GrandBudapest2", type = "continuous"))
    names(pal_COLORS) <-  unique(merge.melt[[tax_level]])
    low_letter <- tolower(substr(tax_level, start = 1, stop = 1))
    pal_COLORS[which(names(pal_COLORS)==paste0(low_letter, "__Unknown"))] <- unknownColor 
    pal_COLORS[which(names(pal_COLORS)=="NA")] <- "green" 
    pal_COLORS[which(names(pal_COLORS)== sprintf("< %d perc. abund.", collapseAbundPerc))] <- collapseColor 
    colScale <- scale_colour_manual(name = tax_level, values = pal_COLORS)
  }
  if (colorPal == "wes3"){
    pal_COLORS<- as.vector(wesanderson::wes_palette(known_col_length, 
                                                    name = "Rushmore1", type = "continuous"))
    names(pal_COLORS) <-  unique(merge.melt[[tax_level]])
    low_letter <- tolower(substr(tax_level, start = 1, stop = 1))
    pal_COLORS[which(names(pal_COLORS)==paste0(low_letter, "__Unknown"))] <- unknownColor 
    pal_COLORS[which(names(pal_COLORS)=="NA")] <- "green" 
    pal_COLORS[which(names(pal_COLORS)== sprintf("< %d perc. abund.", collapseAbundPerc))] <- collapseColor 
    colScale <- scale_colour_manual(name = tax_level, values = pal_COLORS)
  }
  if (colorPal == "vir") { #viridis
    pal_COLORS <- viridis(known_col_length)
    names(pal_COLORS) <-  unique(merge.melt[[tax_level]])
    low_letter <- tolower(substr(tax_level, start = 1, stop = 1))
    pal_COLORS[which(names(pal_COLORS)==paste0(low_letter, "__Unknown"))] <- unknownColor 
    pal_COLORS[which(names(pal_COLORS)=="NA")] <- "green" 
    pal_COLORS[which(names(pal_COLORS)== sprintf("< %d perc. abund.", collapseAbundPerc))] <- collapseColor 
    colScale <- scale_colour_manual(name = tax_level, values = pal_COLORS) # if not using fill
  }
  if (colorPal == "hotcold") { #viridis
    #cold
    pal_COLORS <- viridis(known_col_length)
    reColorList1<- as.vector(wesanderson::wes_palette(Unk_col_length, 
                                                      name = "GrandBudapest1", type = "continuous"))
    # pal_COLORS <- append(pal_COLORS, reColorList1, after = length(pal_COLORS))
    
    names(pal_COLORS) <-  unique(merge.melt[[tax_level]])
    low_letter <- tolower(substr(tax_level, start = 1, stop = 1))
    pal_COLORS[which(names(pal_COLORS)==paste0(low_letter, "__Unknown"))] <- unknownColor 
    pal_COLORS[which(names(pal_COLORS)=="NA")] <- "green" 
    pal_COLORS[which(names(pal_COLORS)== sprintf("< %d perc. abund.", collapseAbundPerc))] <- collapseColor
    colScale <- scale_colour_manual(name = tax_level, values = pal_COLORS) # if not using fill
    #hot
    for (taxName in modified_OTU_names){
      pal_COLORS[which(names(pal_COLORS)==taxName)] <- reColorList1[match(taxName,modified_OTU_names)]
    }
  }
  if (!is.null(reColorTax)){ #if not null
    # for (taxName in reColorTax){ #for each entry
    if (!is.null(reColorList)){ #if list exists
      if (length(reColorList)==length(reColorTax)){ #if enough colors are provided proceed
        for (taxName in reColorTax){
          pal_COLORS[which(names(pal_COLORS)==taxName)] <- reColorList[match(taxName,reColorTax)] 
        }
      }
      else{
        cat("reColorTax and reColorList are not of equal length\n")
        cat("aborting recolor\n")
      }
    }
    else { #reColorList is still NULL
      cat("no colors provided to reColorList, assigning random colors\n")
      require(randomcoloR);
      for (taxName in reColorTax){
        pal_COLORS[which(names(pal_COLORS)==taxName)] <- randomcoloR::randomColor()
      }
    }
  }
  ###
  #draw plot
  gg <- ggplot(merge.melt, 
               aes(x=merge.melt[[group]], y= Abundance,
                   fill=factor(merge.melt[[tax_level]]))) +
    ggtitle(title) + geom_bar(aes(     ),stat = "identity")+ xlab(group) + 
    guides(fill=guide_legend(ncol=Legendcols)) +
    scale_fill_manual(values = pal_COLORS)
  
  if (legend == T){
    cat("plotting with legend\n")
    if (!is.null(facetGrid)){
      cat("plotting with facet_grid ~", facetGrid,"\n")
      gg <- gg + facet_grid(cols = vars(merge.melt[[facetGrid]]), drop = F, scales = "free", space = "free_x") +
        theme(legend.text = element_text(size = 7, margin = margin(r = 0, unit = 'pt')), 
              legend.key=element_blank(),legend.title = element_blank(),
              legend.position=legendLoc,legend.spacing.x = unit(0.1, 'cm'), #distance from legend block
              legend.margin = margin(0,0,0,0),
              axis.text.x = element_text(size = 7, angle=-39, vjust=0.2, hjust = 0.1)) #distance from bars
    }
    else {
      gg <- gg + theme(legend.text = element_text(size = 7, margin = margin(r = 0, unit = 'pt')), 
                       legend.key=element_blank(),legend.title = element_blank(),
                       legend.position=legendLoc,legend.spacing.x = unit(0.1, 'cm'), #distance from legend block
                       legend.margin = margin(0,0,0,0), 
                       axis.text.x = element_text(size = 7, angle=-39, vjust=0.2, hjust = 0.1)) #distance from bars
    }
  } 
  else {
    cat("plotting without legend\n")
    if (!is.null(facetGrid)){
      cat("plotting with facet_grid ~", facetGrid,"\n")
      gg <- gg + facet_grid(cols = vars(merge.melt[[facetGrid]]), drop = F, scales = "free", space = "free_x") +
        theme(legend.position="none", axis.text.x = element_text(size = 7, angle=-39, vjust=0.2, hjust = 0.1))
    }
    else {
      gg <- gg + theme(legend.position="none", axis.text.x = element_text(size = 7, angle=-39, vjust=0.2, hjust = 0.1))
    }
  }
  if (barColSeperator==T){
    #adding black seperators to the plot
    cat("plotting with seperator\n")
    gg <- gg + geom_bar(aes(),stat = "identity", color = "black")
  }
  return(gg)
}



#### Data import and standardization ####

#### plot using microbiome package Genus ####
# from microbiome::plot_taxa_prevalence
check_phyloseq <- function (x, fill_na_taxa = TRUE) {
						# Sanity checks for a phyloseq object. Required with some methods.
						if (!taxa_are_rows(x)) {
							x@otu_table <- otu_table(t(otu_table(x)), taxa_are_rows = TRUE)
						}
						if (fill_na_taxa || is.character(fill_na_taxa)) {
							M <- as.matrix(tax_table(x))
							if (!taxa_are_rows(x)) {
								M <- t(M)
							}
							if (!is.character(fill_na_taxa)) {
								fill_na_taxa <- "Unknown"
							}
							M[is.na(M)] <- fill_na_taxa
							x@tax_table <- tax_table(M)
						}
						x
					}
plot_taxa_prevalence_genus = function (x, level, detection = 0) 
{
	require(microbiome)
	abundance <- NULL
	prevalence <- NULL
	x <-check_phyloseq(x, fill_na_taxa = TRUE)
	if (level == "Phylum") {
		tax.abun <- apply(abundances(x), 1, mean)
		tax.prev <- prevalence(x, detection = detection)
		Phylum <- as.vector(data.frame(tax_table(x))$Phylum)
		Phylum <- as.vector(Phylum)
		Phylum <- as.factor(Phylum)
		xdf <- data.frame(abundance = log(tax.abun), prevalence = tax.prev, 
											Phylum)
		plot.phylum <- ggplot(xdf, aes(x = abundance, y = prevalence, 
																	 color = Phylum)) + geom_point(shape = 16, alpha = 0.9) + 
			xlab("Average count abundance (log scale)") + ylab("Taxa prevalence") + 
			theme_bw() + theme(plot.background = element_blank(), 
												 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
												 panel.background = element_blank(), axis.text.x = element_text(angle = 90, 
												 																															 vjust = 0.5, size = 6)) + geom_vline(xintercept = 0, 
												 																															 																		 linetype = "dashed", color = "grey") + facet_wrap(~Phylum) + 
			theme(strip.background = element_rect(fill = "white"))
		return(plot.phylum)
	}
	else if (level == "Genus") {
		tax.abun <- apply(abundances(x), 1, mean)
		tax.prev <- rowSums(abundances(x) != 0)/nsamples(x)
		Genus <- as.vector(data.frame(tax_table(x))$Genus)
		Genus <- as.vector(Genus)
		Genus <- as.factor(Genus)
		xdf <- data.frame(abundance = log(tax.abun), prevalence = tax.prev, 
											Genus)
		plot.gen <- ggplot(xdf, aes(x = abundance, y = prevalence, 
																color = Genus)) + geom_point(shape = 16, alpha = 0.9) + 
													xlab("Average count abundance (log scale)") + ylab("Taxa prevalence") + 
												 theme_bw() + theme(plot.background = element_blank(), 
												 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
												 panel.background = element_blank(), 
												 axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6)) + 
										geom_vline(xintercept = 0, linetype = "dashed", color = "grey") + 
			facet_wrap(~Genus) + 
			theme(strip.background = element_rect(fill = "white"))
		return(plot.gen)
	}
	else if (level == "Family") {
		tax.abun <- apply(abundances(x), 1, mean)
		tax.prev <- rowSums(abundances(x) != 0)/nsamples(x)
		Family <- as.vector(data.frame(tax_table(x))$Family)
		Family <- as.vector(Family)
		Family <- as.factor(Family)
		xdf <- data.frame(abundance = log(tax.abun), prevalence = tax.prev, 
											Family)
		plot.fam <- ggplot(xdf, aes(x = abundance, y = prevalence, 
																color = Family)) + geom_point(shape = 16, alpha = 0.9) + 
			xlab("Average count abundance (log scale)") + ylab("Taxa prevalence") + 
			theme_bw() + theme(plot.background = element_blank(), 
												 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
												 panel.background = element_blank(), 
												 axis.text.x = element_text(angle = 90,  vjust = 0.5, size = 6)) + 
			geom_vline(xintercept = 0,linetype = "dashed", color = "grey") + 
			facet_wrap(~Family) + 
			theme(strip.background = element_rect(fill = "white"))
		return(plot.fam)
	}
	else if (level == "Order") {
		tax.abun <- apply(abundances(x), 1, mean)
		tax.prev <- rowSums(abundances(x) != 0)/nsamples(x)
		Order <- as.vector(data.frame(tax_table(x))$Order)
		Order <- as.vector(Order)
		Order <- as.factor(Order)
		xdf <- data.frame(abundance = log(tax.abun), prevalence = tax.prev, 
											Order)
		plot.order <- ggplot(xdf, aes(x = abundance, y = prevalence, 
																	color = Order)) + geom_point(shape = 16, alpha = 0.9) + 
			xlab("Average count abundance (log scale)") + ylab("Taxa prevalence") + 
			theme_bw() + theme(plot.background = element_blank(), 
												 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
												 panel.background = element_blank(), axis.text.x = element_text(angle = 90, 
												 																															 vjust = 0.5, size = 6)) + geom_vline(xintercept = 0, 
												 																															 																		 linetype = "dashed", color = "grey") + facet_wrap(~Order) + 
			theme(strip.background = element_rect(fill = "white"))
		return(plot.order)
	}
	else if (level == "Class") {
		tax.abun <- apply(abundances(x), 1, mean)
		tax.prev <- rowSums(abundances(x) != 0)/nsamples(x)
		Class <- as.vector(data.frame(tax_table(x))$Class)
		Class <- as.vector(Class)
		Class <- as.factor(Class)
		xdf <- data.frame(abundance = log(tax.abun), prevalence = tax.prev, 
											Class)
		plot.class <- ggplot(xdf, aes(x = abundance, y = prevalence, 
																	color = Class)) + geom_point(shape = 16, alpha = 0.9) + 
			xlab("Average count abundance (log scale)") + ylab("Taxa prevalence") + 
			theme_bw() + theme(plot.background = element_blank(), 
												 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
												 panel.background = element_blank(), axis.text.x = element_text(angle = 90, 
												 																															 vjust = 0.5, size = 6)) + geom_vline(xintercept = 0, 
												 																															 																		 linetype = "dashed", color = "grey") + facet_wrap(~Class) + 
			theme(strip.background = element_rect(fill = "white"))
		return(plot.class)
	}
}