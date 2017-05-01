## Introduction

Plasmids are mobile DNA elements (MBEs), which allow for rapid evolution and adaption of bacteria to new niches through horizontal transmission of novel traits to different genetic backgrounds. Plasmids and other MBE’s are of great concern to public health since they can encode for traits such as antimicrobial resistance (AMR), toxin production, and virulence. The mobile nature of plasmids can allow for rapid dissemination of AMR through a population, and many plasmids are readily transmissible. Understanding the transmission dynamics of AMR plasmids is critical for public health professionals in order to accurately respond to the rising threat of multi-drug resistant bacteria. Mobility is key to the long-term survival of plasmids and understanding the epidemiology of plasmid-borne traits requires insight into how these vectors are transmitted between hosts. Classification schemes for plasmids based on their mobility have been the subject of several papers which have broadly grouped plasmids based on their relaxase gene and type IV coupling protein (T4CP) into the following categories i) conjugative ii) mobilizable iii) non-mobilizable. Six major relaxase families (MOB) have been previously been identified by relaxase typing based on specific sequence motifs. Currently, there are no tools available to type plasmid assemblies automatically according to MOB groups and predict their mobility. This includes tools such as Plasmid Finder https://cge.cbs.dtu.dk/services/PlasmidFinder/, which focus on the typing of plasmids within draft bacterial genomes, based on the detection of the replication gene(s). Plasmid family information derived from replicon typing can provide information on host range, frequency of AMR elements, among others but it cannot identify whether a given member of a plasmid family possesses the necessary components to transfer itself into another host through conjugation. 

## Release 1.0

Mob-typer will provide in silico predictions of the replicon family, relaxase type, mate-pair formation type and predicted transferability of the plasmid

## Dependancies

Prokka https://github.com/tseemann/prokka
Blast+ https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

## Contact

James Robertson - james.robertson@phac-aspc.gc.ca

## License

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.