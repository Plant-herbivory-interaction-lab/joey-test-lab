library(tidyverse)
library(readxl)
library(lme4)
library(emmeans)
library(car)
library(vegan)
library(devtools)
library(spaa)
library(fastDummies)
library(janitor)
library(DHARMa)
library(glmmTMB)
library(ggeffects)
library(MuMIn)

library(ape)
library(phytools)
library(pez)
library(picante)
library(FD)
library(bipartite)

plants_data_joey<-read_excel("Pollination_neighborhood_data21.xlsx",sheet="plants")
visitors_data<-read_excel("Pollination_neighborhood_data21.xlsx",sheet="visitors")
plots_data<-read_excel("Pollination_neighborhood_data21.xlsx",sheet="plots")
plant_traits<-read_excel("Pollination_neighborhood_data21.xlsx",sheet="plant_spp_traits")
visitor_traits<-read_excel("Pollination_neighborhood_data21.xlsx",sheet="visitor_spp_traits")

###subset and process plants data
plants<-full_join(plants_data,plant_traits,by="flower_sp")
plants_v<-plants%>%select(site,date,plot,flower_no,flower_sp)%>%
  group_by(site,date,plot,flower_sp)%>%
  pivot_wider(names_from = flower_sp, values_from=flower_no,values_fn=mean,values_fill=0)
plants_v$total<-rowSums(plants_v[4:29])#check decimal
plants_v$monfreq<-plants_v$`Monarda fistulosa`/plants_v$total
plants_v$foc_dens<-plants_v$`Monarda fistulosa`
sppmat<-as.matrix(plants_v[4:29])
plotmat<-cbind(plants_v[1:3],plants_v[31:32])
plotmat$shannon<-diversity(sppmat)
plotmat$richness<-specnumber(sppmat)
plotmat$evenness<-plotmat$shannon/log(plotmat$richness)
sppmat_n<-as.matrix(plants_v[5:29])
plotmat$shannon_n<-diversity(sppmat_n)
plotmat$richness_n<-specnumber(sppmat_n)
plotmat$evenness_n<-plotmat$shannon_n/log(plotmat$richness_n)
plot_full<-full_join(plotmat,plots_data)#add monfreq to rest

###subset and process insects data
insects<-full_join(visitors_data,visitor_traits,by="visitor_sp")
insects$visitor_sp<-gsub(" ","_",insects$visitor_sp)
insects_v<-insects%>%filter(flower_sp=="Monarda fistulosa")%>%
  select(site,date,plot,visitor_no,visitor_sp)%>%
  group_by(site,date,plot,visitor_sp)%>%
  pivot_wider(names_from = visitor_sp, values_from=visitor_no,values_fn=mean,values_fill=0)
sppmat2<-as.matrix(insects_v[4:27])
plotmat2<-insects_v[1:3]
plotmat2$shannon_vis<-diversity(sppmat2)
plotmat2$richness_vis<-specnumber(sppmat2)
plotmat2$evenness_vis<-plotmat2$shannon_vis/log(plotmat2$richness_vis)

###neighbor height
nheight<-plants_data%>%select(site,date,plot,flower_sp,flower_tall)%>%
  group_by(site,date,plot,flower_sp)%>%
  pivot_wider(names_from = flower_sp, values_from=flower_tall,values_fn=mean)
nheight$nmax<-apply(nheight[5:29],1,max,na.rm=T)
nheight$nmean<-apply(nheight[5:29],1,mean,na.rm=T)
nheight[mapply(is.infinite,nheight)]<-0
nheight[mapply(is.nan,nheight)]<-0
nheight$tallest<-ifelse(nheight$nmax>nheight$`Monarda fistulosa`,"neighbor","focal")
nheight$tallest2<-ifelse(nheight$nmean>nheight$`Monarda fistulosa`,"neighbor","focal")
nheight2<-nheight%>%select(site,date,plot,tallest,tallest2,`Monarda fistulosa`)
nheight[mapply(is.infinite,nheight)]<-0

plot_full<-full_join(plotmat,plots_data)
plot_full_a<-full_join(plot_full,plotmat2)
plot_full2<-full_join(nheight2,plot_full_a)

plants_insects<-full_join(plants,insects,by=c("site","date","plot","flower_sp"))
full_data<-full_join(plants_insects,plots_data,by=c("site","date","plot"))

#neighbor data
ndat1<-full_data%>%filter(flower_sp!="Monarda fistulosa")%>%
  group_by(site,date,plot)%>%
  summarise(prop_yellow=mean(yellow,na.rm=T),num_yellow=sum(yellow,na.rm=T))%>%
  mutate(has_yellow=ifelse(prop_yellow>0,"yellow_neighbor","no_yellow_neighbor"))
ndat2<-full_data%>%filter(flower_sp!="Monarda fistulosa")%>%
  group_by(site,date,plot)%>%
  summarise(prop_pink=mean(pink,na.rm=T),num_pink=sum(pink,na.rm=T),nbr_dens=sum(flower_no,na.rm=T))%>%
  mutate(has_pink=ifelse(prop_pink>0,"pink_neighbor","no_pink_neighbor"))
ndat<-full_join(ndat1,ndat2)

mondat<-full_data%>%filter(flower_sp=="Monarda fistulosa")%>%
  full_join(.,plot_full2)
mondat2a<-mondat%>%
  group_by(site,date,plot,monfreq,foc_dens,shannon,richness,evenness,shannon_n,richness_n,evenness_n,
           shannon_vis,richness_vis,evenness_vis,flower_no,tallest,tallest2)%>%
  summarise(visitor_no=sum(visitor_no,na.rm = T))
mondat2<-full_join(mondat2a,ndat)
mondat2$percap_vis<-mondat2$visitor_no/mondat2$flower_no

###plot flower diversity vs plot monfis visitor diversity
ggplot(data=mondat2,mapping=aes(x=shannon,y=shannon_vis))+
  geom_jitter(aes(color=flower_no))+
  geom_smooth(method="lm",color="black")+
  labs(x="floral diversity (Shannon index)",y="MonFis visitor diversity")+
  theme_bw()+
  scale_color_gradient2(high="purple",mid="gray",low="green",midpoint=30)
median(mondat2$flower_no)
mean(mondat2$flower_no)

###plot monfis flower number vs plot monfis visitor diversity
ggplot(data=mondat2,mapping=aes(x=flower_no,y=shannon_vis))+
  geom_jitter()+
  geom_smooth(method="lm",color="black")+
  labs(x="MonFis inflorescence no.",y="MonFis visitor diversity")+
  theme_bw()

###plot flower diversity vs plot per capita monfis visitors
ggplot(data=mondat2,mapping=aes(x=shannon,y=percap_vis))+
  geom_point(size=3,alpha=0.5,show.legend=F)+
  geom_smooth(method="lm",color="black")+
  labs(x="Floral diversity (Shannon index)",y="Visitors \n per Monfis inflorescence")+
  theme_bw(base_size=15)

vis1<-glmer(percap_vis~shannon+(1|site/plot),data=mondat2,family=gaussian())#result H1a (all)
Anova(vis1)
summary(vis1)
r.squaredGLMM(vis1)
emtrends(vis1,~shannon,var="shannon",type="response")
vis_r<-glmer(percap_vis~richness+(1|site/plot),data=mondat2,family=gaussian())#result H1b (all)
Anova(vis_r)
summary(vis_r)
r.squaredGLMM(vis_r)
emtrends(vis_r,~richness,var="richness",type="response")
vis2<-glmer(percap_vis~shannon+(1|site),data=mondat2[which(mondat2$percap_vis>0),],family=gaussian())
Anova(vis2)

hist(mondat2$percap_vis)
hist(log(mondat2$percap_vis))

ggplot(data=mondat2,mapping=aes(x=nbr_dens,y=percap_vis))+
  geom_point(size=3,aes(color=monfreq),show.legend=F)+
  geom_smooth(color="black",method="lm")+
  labs(x="neighboring floral units per quadrat",y="visitors \n per Monfis inflorescence")+
  theme_bw(base_size=15)+
  scale_color_gradient2(high="purple",mid="gray",low="green",midpoint=.5,name="MonFis \n frequency")
hist(mondat2$nbr_dens)
hist(log(mondat2$nbr_dens))
vis3<-lmer(percap_vis~nbr_dens+(1|site/plot),data=mondat2)#result H3(all)
Anova(vis3)
summary(vis3)
r.squaredGLMM(vis3)
emtrends(vis3,~nbr_dens,var="nbr_dens",type="response")
plot(simulateResiduals(vis3))
vis3b<-lmer(percap_bom~nbr_dens+(1|site/plot),data=mondat_b2)#result H3(bombus)
Anova(vis3b)
summary(vis3b)
r.squaredGLMM(vis3b)
emtrends(vis3b,~nbr_dens,var="nbr_dens",type="response")
vis3c<-lmer(percap_bom~foc_dens+(1|site/plot),data=mondat_b2)
Anova(vis3c)
summary(vis3c)
r.squaredGLMM(vis3c)
emtrends(vis3c,~foc_dens,var="foc_dens",type="response")
vis3d<-lmer(percap_vis~foc_dens+(1|site/plot),data=mondat2)
Anova(vis3d)
summary(vis3d)
r.squaredGLMM(vis3d)
emtrends(vis3d,~foc_dens,var="foc_dens",type="response")

###plot monfis flower number vs plot per capita monfis visitors
ggplot(data=mondat2,mapping=aes(x=flower_no,y=percap_vis))+
  geom_jitter()+
  geom_smooth(method="lm",color="black")+
  labs(x="MonFis inflorescence no.",y="MonFis visitors per inflorescence")+
  theme_bw()

###subset and process bombus data
bombus_v<-insects%>%filter(flower_sp=="Monarda fistulosa",visitor_gen=="Bombus")%>%
  select(site,date,plot,visitor_no,visitor_sp)%>%
  group_by(site,date,plot,visitor_sp)%>%
  pivot_wider(names_from = visitor_sp, values_from=visitor_no,values_fn=mean,values_fill=0)
sppmat3<-as.matrix(bombus_v[4:10])
plotmat3<-bombus_v[1:3]
plotmat3$shannon_bom<-diversity(sppmat3)
plotmat3$richness_bom<-specnumber(sppmat3)
plotmat3$evenness_bom<-plotmat3$shannon_bom/log(plotmat3$richness_bom)

plot_full_b<-full_join(plot_full,plotmat3)
plot_full3<-full_join(nheight2,plot_full_b)

mondat_b<-full_data%>%filter(flower_sp=="Monarda fistulosa",visitor_gen=="Bombus")%>%
  full_join(.,plot_full3)
mondat_b2a<-mondat_b%>%
  group_by(site,date,plot,monfreq,foc_dens,shannon,richness,evenness,shannon_n,richness_n,
           evenness_n,shannon_bom,richness_bom,evenness_bom,flower_no,tallest,tallest2)%>%
  summarise(visitor_no=sum(visitor_no,na.rm = T))
mondat_b2<-full_join(mondat_b2a,ndat)
mondat_b2$percap_bom<-mondat_b2$visitor_no/mondat_b2$flower_no
mondat_b2$percap_bom[is.na(mondat_b2$percap_bom)]<-0
mondat_b2$has_pink[is.na(mondat_b2$has_pink)]<-"no neighbor"
mondat2$has_pink[is.na(mondat2$has_pink)]<-"no neighbor"
mondat_b2$has_yellow[is.na(mondat_b2$has_yellow)]<-"no neighbor"
mondat2$has_yellow[is.na(mondat2$has_yellow)]<-"no neighbor"
mondat_b2$pinkfreq<-mondat_b2$num_pink/mondat_b2$nbr_dens
mondat_b2$yellowfreq<-mondat_b2$num_yellow/mondat_b2$nbr_dens

ggplot(data = mondat_b2,mapping=aes(x=has_pink,y=percap_bom))+
  geom_boxplot(outlier.shape=NA,size=.8)+
  geom_jitter(width=.2,alpha=0.5,size=3,show.legend=F)+
  labs(x="Neighbor color",y="Bombus spp. visitors \n per MonFis inflorescence")+
  scale_x_discrete(labels=c("no neighbor"="no neighbor","no_pink_neighbor"="no pink-purple \n neighbor","pink_neighbor"="pink-purple \n neighbor"))+
  theme_bw(base_size=20)
#fig3
pinkmod<-lmer(percap_bom~has_pink+(1|site/plot),data=mondat_b2)
Anova(pinkmod)#result H2(bombus)
plot(simulateResiduals(pinkmod))
emmeans(pinkmod,pairwise~has_pink)
r.squaredGLMM(pinkmod)
pinkmod_all<-lmer(percap_vis~has_pink+(1|site/plot),data=mondat2)
Anova(pinkmod_all)#result H2(all)
plot(simulateResiduals(pinkmod_all))
emmeans(pinkmod_all,pairwise~has_pink)
r.squaredGLMM(pinkmod_all)
yellowmod<-lmer(percap_bom~has_yellow+(1|site),data=mondat_b2)
Anova(yellowmod)
yellowmod_all<-lmer(percap_vis~has_yellow+(1|site),data=mondat2)
Anova(yellowmod_all)


ggplot(data=mondat_b2,mapping=aes(x=prop_pink,y=percap_bom))+
  geom_point(size=3,color="purple")+
  geom_smooth(color="black",method="lm",formula=y~x+I(x^2))+
  labs(x="proportion of pink/purple flowers in quadrat",y="bumblebee visitors \n per MonFis inflorescence")+
  theme_bw(base_size=15)

ggplot(data=mondat_b2,mapping=aes(x=prop_yellow,y=percap_bom))+
  geom_point(size=3,color="goldenrod")+
  geom_smooth(color="black",method="lm",formula=y~x+I(x^2))+
  labs(x="proportion of yellow flowers in quadrat",y="bumblebee visitors \n per MonFis inflorescence")+
  theme_bw(base_size=15)

ggplot(data=mondat_b2,mapping=aes(x=num_pink,y=prop_pink))+
  geom_point(size=3,color="purple")+
  geom_smooth(color="black",method="lm")+
  labs(x="neighboring pink/purple floral units per quadrat",y="proportion of neighbors w/ pink/purple flowers")+
  theme_bw(base_size=15)

ggplot(data=mondat_b2,mapping=aes(x=num_yellow,y=prop_yellow))+
  geom_point(size=3,color="goldenrod")+
  geom_smooth(color="black",method="lm")+
  labs(x="neighboring yellow floral units per quadrat",y="proportion of neighbors w/ pink/purple flowers")+
  theme_bw(base_size=15)

###plot flower diversity vs plot per capita monfis visitors (bombus)
ggplot(data=mondat_b2,mapping=aes(x=shannon,y=percap_bom))+
  geom_point(size=3,alpha=0.5,show.legend=F)+
  geom_smooth(method="lm",color="black")+
  labs(x="Floral diversity (Shannon index)",y="Bombus spp. visitors \n per MonFis inflorescence")+
  theme_bw(base_size=20)
#fig1

hist(mondat_b2$percap_bom)
hist(mondat_b2$visitor_no)
bomvis<-lmer(percap_bom~shannon+(1|site),data=mondat_b2)
Anova(bomvis)
plot(simulateResiduals(bomvis))
bomvis1<-lmer(percap_bom~shannon+(1|site/plot),data=mondat_b2)#result H1a(Bombus)
Anova(bomvis1)
plot(simulateResiduals(bomvis1))
summary(bomvis1)
r.squaredGLMM(bomvis1)
emtrends(bomvis1,~shannon,var="shannon",type="response")

bomvis2<-glmmTMB(visitor_no~shannon+(1|site/plot/date),data=mondat_b2,family=poisson,offset=log(flower_no))
Anova(bomvis2)
plot(simulateResiduals(bomvis2))
summary(bomvis2)
bomvis3<-glmer(percap_bom~shannon+(1|site),data=mondat_b2[which(mondat_b2$percap_bom>0),],family=gaussian())
Anova(bomvis3)

ggplot(data=mondat_b2,mapping=aes(x=shannon_n,y=percap_bom))+
  geom_jitter(size=3,aes(color=monfreq),show.legend=T)+
  geom_smooth(method="lm",color="black")+
  labs(x="neighbor floral diversity (Shannon index)",y="bumblebee visitors \n per beebalm inflorescence")+
  theme_bw(base_size=20)+
  scale_color_gradient2(high="purple",mid="grey",low="green",midpoint=.5,name="beebalm \n frequency")

bomvis4<-lmer(percap_bom~shannon_n+(1|site/plot),data=mondat_b2)
Anova(bomvis4)
plot(simulateResiduals(bomvis4))
bomvis4a<-lmer(log(percap_bom+1)~shannon_n+(1|site/plot/date),data=mondat_b2)
Anova(bomvis4a)
plot(simulateResiduals(bomvis4a))
bomvis4b<-glmmTMB(visitor_no~shannon_n+(1|site/plot/date),data=mondat_b2,family="poisson")
Anova(bomvis4b)
plot(simulateResiduals(bomvis4b))
bomvis4c<-glmmTMB(visitor_no~shannon_n*monfreq+(1|site/plot/date),data=mondat_b2,family="poisson")
Anova(bomvis4c)
plot(simulateResiduals(bomvis4c))

ggplot(data=mondat_b2[which(mondat_b2$percap_bom>0),],mapping=aes(x=shannon,y=percap_bom))+
  geom_jitter(size=3,aes(),show.legend=T)+
  geom_smooth(method="lm",color="black")+
  facet_wrap(vars(site))

ggplot(data=mondat_b2,mapping=aes(x=richness,y=percap_bom))+
  geom_point(size=3,alpha=0.5,show.legend=F)+
  geom_smooth(method="lm",color="black")+
  labs(x="Floral richness",y="Bombus spp. visitors \n per MonFis inflorescence")+
  theme_bw(base_size=20)
#fig2

bomvis_r<-lmer(percap_bom~richness+(1|site/plot),data=mondat_b2)#result H1b(Bombus)
Anova(bomvis_r)
plot(simulateResiduals(bomvis_r))
emtrends(bomvis_r,~richness,var="richness",type="response")
r.squaredGLMM(bomvis_r)
bomvis_r1<-lmer(percap_bom~richness+(1|site/plot/date),data=mondat_b2)
Anova(bomvis_r1)
plot(simulateResiduals(bomvis_r1))
bomvis_r2<-glmmTMB(visitor_no~richness+(1|site/plot/date),data=mondat_b2,family="poisson")
Anova(bomvis_r2)
plot(simulateResiduals(bomvis_r2))
bomvis_r3<-glmmTMB(visitor_no~richness*foc_dens+(1|site/plot/date),data=mondat_b2,family="poisson")
Anova(bomvis_r3)
plot(simulateResiduals(bomvis_r3))

bomvis_r3<-lmer(percap_bom~richness_n+(1|site/plot),data=mondat_b2)
Anova(bomvis_r3)
plot(simulateResiduals(bomvis_r3))
bomvis_r3a<-glmmTMB(visitor_no~richness_n+(1|site/plot),data=mondat_b2,family="poisson")
Anova(bomvis_r3a)
plot(simulateResiduals(bomvis_r3a))
bomvis_r3b<-glmmTMB(visitor_no~richness_n*monfreq+(1|site/plot),data=mondat_b2,family="poisson")
Anova(bomvis_r3b)
plot(simulateResiduals(bomvis_r3b))
bomvis_r3b_pred<-ggpredict(bomvis_r3b,terms="richness_n")
plot(bomvis_r3b_pred,add.data=T)

###plot monfis flower number vs plot per capita monfis visitors (bombus)
ggplot(data=mondat_b2,mapping=aes(x=flower_no,y=visitor_no))+
  geom_point(color="purple",size=3,show.legend=F)+
  geom_smooth(method="lm",formula=y~x+I(x^2),color="black")+
  #scale_x_continuous(trans="log2")+
  labs(x="beebalm inflorescence no.",y="bumblebee visitors to beebalm")+
  theme_bw(base_size=20)

ggplot(data=mondat_duf,mapping=aes(x=flower_no,y=visitor_no))+
  geom_point(color="purple",size=3,show.legend=F)+
  geom_smooth(method="lm",formula=y~x+I(x^2),color="black")+
  #scale_x_continuous(trans="log2")+
  labs(x="beebalm inflorescence no.",y="D. monardae visitors to beebalm")+
  theme_bw(base_size=20)

#Halictidae (pollen foraging only)
mondat_hal<-full_data%>%filter(flower_sp=="Monarda fistulosa",visitor_fam=="Halictidae")%>%
  full_join(.,plot_full2)
mondat_hal<-mondat_hal%>%
  group_by(site,date,plot,monfreq,foc_dens,shannon,richness,evenness,
           shannon_n,richness_n,evenness_n,flower_no,tallest,tallest2)%>%
  summarise(visitor_no=sum(visitor_no,na.rm = T))
mondat_hal<-full_join(mondat_hal,ndat)
mondat_hal$percap_hal<-mondat_hal$visitor_no/mondat_hal$flower_no
mondat_hal$percap_hal[is.na(mondat_hal$percap_hal)]<-0

ggplot(data=mondat_hal,mapping=aes(x=shannon_n,y=percap_hal))+
  geom_point(size=3,aes(color=monfreq),show.legend=T)+
  geom_smooth(method="lm",color="black")+
  labs(x="floral diversity (Shannon index)",y="sweat bee visitors \n per MonFis inflorescence")+
  theme_bw(base_size=15)+
  scale_color_gradient2(high="purple",mid="grey",low="green",midpoint=.5,name="MonFis \n frequency")

halvis<-glmer(percap_hal~shannon+(1|site),data=mondat_hal,family=gaussian())
Anova(halvis)

#Dufourea monardae (specialist)
mondat_duf<-full_data%>%filter(flower_sp=="Monarda fistulosa",visitor_gen=="Dufourea")%>%
  full_join(.,plot_full2)
mondat_duf<-mondat_duf%>%
  group_by(site,date,plot,monfreq,foc_dens,shannon,richness,evenness,
           shannon_n,richness_n,evenness_n,flower_no,tallest,tallest2)%>%
  summarise(visitor_no=sum(visitor_no,na.rm = T))
mondat_duf<-full_join(mondat_duf,ndat)
mondat_duf$percap_duf<-mondat_duf$visitor_no/mondat_duf$flower_no
mondat_duf$percap_duf[is.na(mondat_duf$percap_duf)]<-0

ggplot(data=mondat_duf,mapping=aes(x=shannon,y=percap_duf))+
  geom_jitter(size=3,aes(color=monfreq),show.legend=T)+
  geom_smooth(method="lm",color="black")+
  labs(x="floral diversity (Shannon index)",y="D. monardae \n per MonFis inflorescence")+
  theme_bw(base_size=15)+
  scale_color_gradient2(high="purple",mid="grey",low="green",midpoint=.5,name="MonFis \n frequency")

dufvis<-glmer(percap_duf~shannon+(1|site),data=mondat_duf,family=gaussian())
Anova(dufvis)

###visitor phylogenetic analysis
bees<-insects%>%filter(flower_sp=="Monarda fistulosa",visitor_order=="Hymenoptera")%>%
  filter(visitor_fam!="Crabronidae",visitor_fam!="Sphecidae")

beetree_full<-read.nexus("bee_spp_phyl.nex")
bee_tax<-insects%>%filter(flower_sp=="Monarda fistulosa",visitor_order=="Hymenoptera")%>%
  filter(visitor_fam!="Crabronidae",visitor_fam!="Sphecidae")%>%
  select(visitor_sp,visitor_subgen,visitor_gen,visitor_fam)%>%
  group_by(visitor_sp,visitor_subgen,visitor_gen,visitor_fam)%>%
  filter(visitor_sp!="Augochlorini_sp")
bee_tax$visitor_sp<-gsub(" ","_",bee_tax$visitor_sp)
beegentree<-congeneric.merge(beetree_full,bee_tax$visitor_sp)
beegentree_trim<-keep.tip(beegentree,bee_tax$visitor_sp)
beegentree_trim<-force.ultrametric(beegentree_trim)
plot.phylo(beegentree_trim,type="phylogram",cex=.75,main="Bee spp phylogeny")
beegentree_trim$tip.label

bees_v<-insects%>%filter(flower_sp=="Monarda fistulosa",visitor_order=="Hymenoptera",visitor_fam!="Crabronidae",visitor_fam!="Sphecidae")%>%
  filter(visitor_sp!="Augochlorini")%>%
  select(site,date,plot,visitor_no,visitor_sp)%>%
  group_by(site,date,plot,visitor_sp)%>%
  pivot_wider(names_from = visitor_sp, values_from=visitor_no,values_fn=mean,values_fill=0)
bees_p<-subset(bees_v,select=-c(site,date,plot))
bees_p<-subset(bees_p,select=-Augochlorini_sp)
bees_p<-bees_p%>%select(order(colnames(.))) %>% as.matrix()

mondat_bees<-subset(bees_v,select=c(site,date,plot))
mondat_bees1<-left_join(mondat_bees,mondat2a)

mondat_bees1$pMPD<-mpd(bees_p,cophenetic(beegentree_trim),abundance.weighted=T)
mondat_bees1$pPD<-pd(bees_p,beegentree_trim,include.root=F)[,1]
mondat_bees1$prich<-pd(bees_p,beegentree_trim,include.root=F)[,2]
mondat_bees1$phyrich<-psr(bees_p,beegentree_trim,compute.var=TRUE,scale.vcv=TRUE)


interact_v<-full_data%>%select(flower_sp,visitor_no,visitor_sp)%>%
  group_by(flower_sp,visitor_sp)%>%
  filter(flower_sp!='Dianthus armeria')%>%
  pivot_wider(names_from = flower_sp,values_from=visitor_no,values_fn=sum,values_fill=0)
intmat<-as.matrix(interact_v[2:26])
intmat<-intmat[-8,]
morisita<-niche.overlap(intmat,method="morisita")
int_mormat<-as.matrix(morisita)
int_mormat[is.nan(int_mormat)]<-0
intmat2<-t(intmat)
int_euc<-vegdist(intmat2,method="euclidean",na.rm=T)
int_eucmat<-as.matrix(int_euc)

ptraits_v<-plant_traits%>%select(flower_sp,petal_color,corolla_depth,orientation,flor_shape,symmetry)
ptraits_v<-dummy_cols(ptraits_v,select_columns=c('petal_color','orientation','flor_shape','symmetry'),remove_selected_columns=T)
ptraits_v<-as.data.frame(t(ptraits_v))%>%row_to_names(row_number=1)
ptrtmat<-as.matrix(ptraits_v[1:25])
row.names(ptrtmat)<-NULL
class(ptrtmat)<-"numeric"
ptrt_mori<-niche.overlap(ptrtmat,method="morisita")
ptrt_mormat<-as.matrix(ptrt_mori)
ptrtmat2<-t(ptrtmat)
ptrt_euc<-vegdist(ptrtmat2,method="euclidean")
ptrt_eucmat<-as.matrix(ptrt_euc)

mantel(int_mormat,ptrt_mormat)
mantel(int_eucmat,ptrt_eucmat)

interact_v2<-interact_v[-8,]
interact_v2<-column_to_rownames(interact_v2,var="visitor_sp")
intmat3<-t(as.matrix(interact_v2))
plotweb(intmat3)
visweb(intmat3)
visweb(intmat3,type="diagonal")
networklevel(intmat3)
