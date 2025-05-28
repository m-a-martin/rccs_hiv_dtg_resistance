suppressMessages(library(tidyverse))
suppressMessages(library(tidytree))
suppressMessages(library(treeio))
suppressMessages(require(ggtree))
suppressMessages(require(ape))
suppressMessages(library(shadowtext))
suppressMessages(library(patchwork))
suppressMessages(source('scripts/utils.R'))


muts = bind_rows(
	read_tsv('data/rakai_drug_resistance_mut_all.tsv', show_col_types=FALSE),
	read_tsv('data/other_rakai_drug_resistance_mut.tsv', show_col_types=FALSE))

ins153y = muts %>% 
	filter(mut == 'inS153Y') %>% 
	select(study_id, int_date, mut) %>% 
	unique() %>%
	mutate(ins153y=TRUE) %>%
	select(-mut)

ins153y_co = muts %>% inner_join(ins153y, by=c('study_id', 'int_date')) %>%
	select(study_id, int_date, mut, ins153y) %>% unique() %>%
	mutate(z = TRUE) %>%
	pivot_wider(names_from=mut, values_from=z, values_fill=FALSE) %>%
	pivot_longer(-c(study_id, int_date, ins153y)) %>%
	left_join(read_tsv('config/mut_class.tsv', show_col_types=FALSE), by=c('name'='mut'))

ins153y_co = ins153y_co %>%
	mutate(name = ordered(name, levels=
		(ins153y_co %>% select(class, name) %>% unique() %>%
			mutate(pos = as.numeric(substr(name, 4, nchar(name)-1))) %>%
			arrange(name != 'inS153Y', class, pos))$name))

p = ggplot(ins153y_co, aes(y=study_id, x=name, fill=class, alpha=value)) +
	geom_tile(color='#333333') +
	scale_fill_manual(values=colors, guide=NULL) +
	scale_x_discrete(expand=c(0,0), name='mutation') +
	scale_y_discrete(expand=c(0,0), name='participant', labels=NULL) +
	scale_alpha_manual(values=c(`TRUE`=1,`FALSE`=0), guide=NULL) +
	gtheme +
	theme(panel.grid.major = element_blank(),
		panel.grid.minor = element_blank())


ggsave('figures/pdf/ins153y_co.pdf', p, width=6.4, height=4.8)
ggsave('figures/png/ins153y_co.png', p, width=6.4, height=4.8)
