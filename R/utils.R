
# split gene or neuron list into vector
split_text_to_vector <- function(text){
  stringr::str_split(text,
                     pattern = "[[:space:],;]")[[1]] %>%
    stringi::stri_remove_empty()
}

# ensure we have a list of valid neurons
validate_neurons <- function(neurs, neurons_table){
  
  # replace synonyms
  neurs <- recode(neurs,
                  SENS="SENSORY",
                  INTER="INTERNEURON",
                  INTERNEURONS="INTERNEURON",
                  MOTORNEURON="MOTOR",
                  MOTORNEURONS="MOTOR",
                  MOTONEURON="MOTOR",
                  MOTONEURONS="MOTOR",
                  PHARYNX="PHARYNGEAL",
                  CILIA="CILIATED",
                  CHOLINERGIC="ACH",
                  ACHERGIC="ACH",
                  ACETYLCHOLINERGIC="ACH",
                  ACETYLCHOLINE="ACH",
                  CHO="ACH",
                  SEROTONERGIC="SEROTONIN",
                  SEROTONINERGIC="SEROTONIN",
                  `5-HT`="SEROTONIN",
                  `5HT`="SEROTONIN",
                  GLUTAMATERGIC="GLUTAMATE",
                  GLU="GLUTAMATE",
                  GLUERGIC="GLUTAMATE",
                  GLUTA="GLUTAMATE",
                  GABAERGIC="GABA",
                  OCTO="OCTOPAMINE",
                  OCTOPAMINERGIC="OCTOPAMINE")
  
  
  if(length(neurs) == 1L && stringr::str_to_upper(neurs) == "ALL"){
    return(neurons_table$Neuron_type)
  }
  
  
  # Subcategories
  if(any(neurs == "SENSORY")){
    neurs <- setdiff(neurs, "SENSORY") %>%
      c(neurons_table$Neuron_type[stringr::str_detect(neurons_table$Modality,
                                                      "sensory")])
  }
  if(any(neurs == "INTERNEURON")){
    neurs <- setdiff(neurs, "INTERNEURON") %>%
      c(neurons_table$Neuron_type[stringr::str_detect(neurons_table$Modality,
                                                      "interneuron")])
  }
  if(any(neurs == "MOTOR")){
    neurs <- setdiff(neurs, "MOTOR") %>%
      c(neurons_table$Neuron_type[stringr::str_detect(neurons_table$Modality,
                                                      "motor")])
  }
  if(any(neurs == "PHARYNGEAL")){
    neurs <- setdiff(neurs, "PHARYNGEAL") %>%
      c(neurons_table$Neuron_type[stringr::str_detect(neurons_table$Modality,
                                                      "pharyngeal")])
  }
  if(any(neurs == "CILIATED")){
    neurs <- setdiff(neurs, "CILIATED") %>%
      c(neurons_table$Neuron_type[neurons_table$Ciliated == "yes"])
  }
  if(any(neurs == "ACH")){
    neurs <- setdiff(neurs, "ACH") %>%
      c(neurons_table$Neuron_type[stringr::str_detect(neurons_table$Neurotransmitter,
                                                      "ach")])
  }
  if(any(neurs == "SEROTONIN")){
    neurs <- setdiff(neurs, "SEROTONIN") %>%
      c(neurons_table$Neuron_type[stringr::str_detect(neurons_table$Neurotransmitter,
                                                      "serotonin")])
  }
  if(any(neurs == "GLUTAMATE")){
    neurs <- setdiff(neurs, "GLUTAMATE") %>%
      c(neurons_table$Neuron_type[stringr::str_detect(neurons_table$Neurotransmitter,
                                                      "glutamate")])
  }
  if(any(neurs == "GABA")){
    neurs <- setdiff(neurs, "GABA") %>%
      c(neurons_table$Neuron_type[stringr::str_detect(neurons_table$Neurotransmitter,
                                                      "gaba")])
  }
  if(any(neurs == "OCTOPAMINE")){
    neurs <- setdiff(neurs, "OCTOPAMINE") %>%
      c(neurons_table$Neuron_type[stringr::str_detect(neurons_table$Neurotransmitter,
                                                      "octopamine")])
  }
  
  
  
  
  
  
  # Remove unknown neurons
  unknown_neurs <- setdiff(neurs, neurons_table$Neuron_type)
  
  if(length(unknown_neurs) != 0L){
    warning("Neuron name not recognized: ", unknown_neurs)
  }
  
  setdiff(neurs, unknown_neurs)
}
