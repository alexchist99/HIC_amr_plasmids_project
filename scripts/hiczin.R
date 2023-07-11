#Normalize metagenomic Hi-C data and detect spurious contacts using zero-inflated Negative Binominal regression frameworks
#Auther and maintainer: Yuxuan Du <yuxuandu@usc.edu>
#HiCzin R script depends on 'glmmTMB' package and 'optparse' package

#setwd('~/hic_metagenomics/Reanimation_gut/IC6')
library("optparse")
library("glmmTMB")
library('data.table')
options(scipen = 999)

option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NA,
              help = "input raw metagenomic Hi-C contacts"),
  make_option(c("-s", "--sample"), type = "character", default = NA,
              help = "input samples of the intra-species contacts"),
  make_option(c("-o", "--output"), type = "character", default = NA,
              help = "output normalized Hi-C contacts by HiCzin"),
  make_option(c("-d", "--discard"), type = "character", default = NA,
              help = "output normalized Hi-C contacts that are regarded as spurious contacts and discarded"),
  make_option(c("-f", "--info"), type = "character", default = NA,
              help = "contig information file"),
  make_option(c("-c", "--coverage"), type = "character", default = NA,
              help = "contigs' converage file"),
  make_option(c("-t", "--threshold"), type = "character", default = 0.1,
              help = "preselected threshold for discarding spurious contacts"),
  make_option(c("-u", "--unlabeled"), action = "store_true", default = FALSE,
              help = "use ulabeled mode of HiCzin"),
  make_option(c("-n", "--nb"), action = "store_true", default = FALSE,
              help = "use negative binomial regression"),
  make_option(c("-m", "--model"), type = "character", default = FALSE,
              help = "one of 'unl', 'nb', 'none'")
)


opt_parser <- OptionParser(option_list=option_list)
#print(opt_parser)

opt <- parse_args(opt_parser)
#print(opt)

allcontact_file_name = opt$input
sample_file_name = opt$sample
output_file_name = opt$output
discarded_contact_name = opt$discard
contig_info_name = opt$info
model = opt$model

#print(allcontact_file_name)


# model='unl'
# allcontact_file_name = "links_count_strip.txt"
# sample_file_name = "links_good_strip.txt"
# output_file_name = paste0("good_IC6_",model,".csv")
# contig_info_name = "contig_info.txt"
# #discarded_contact_name = "spur.csv"
# #contig_info_name = "len.txt"
# #coverage_name = "depth.tsv"
thres = 0.1

all_contacts = fread(allcontact_file_name,sep=' ',header = T)
all_contacts = as.data.frame(all_contacts)
colnames(all_contacts) = c('contacts','index1' , 'index2' )

#print(nrow(all_contacts))


sample_data = read.csv(sample_file_name , header = T , sep = ' ')
#sample_data = fread(sample_file_name,sep=' ',header = T)
sample_data = as.data.frame(sample_data)
colnames(sample_data) = c('contacts','index1' , 'index2')
print(head(sample_data))


contig_info = read.csv(contig_info_name , header = F , sep = '\t' )
contig_info = as.data.frame(contig_info)
colnames(contig_info) = c('contig_name','length','coverage')
print(head(contig_info))

sample_len = rep(0 , nrow(sample_data))
#print(sample_len)
sample_cov = rep(0 , nrow(sample_data))
message('processing intra-sample data')
for(i in 1:nrow(sample_data))
{

  sample_len[i] = log(as.numeric(contig_info[contig_info$contig_name==sample_data[i , 'index1'] , 'length']) * 
                        as.numeric(contig_info[contig_info$contig_name==sample_data[i , 'index2'] , 'length']))
  #print(sample_len[i])

 
  sample_cov[i] = log(as.numeric(contig_info[contig_info$contig_name==sample_data[i , 'index1'] , 'coverage']) * 
                        as.numeric(contig_info[contig_info$contig_name==sample_data[i , 'index2'] , 'coverage']))
  #print(sample_cov[i])
}

sampleCon = as.numeric(sample_data[ , 'contacts'])
#print(sampleCon)


#sample_cov = sample_cov[!is.infinite(sample_cov)]
#sample_cov = sample_cov[!is.nan(sample_cov)]
#print(length(which(is.nan(sample_cov))))
print(length(which(is.infinite(sample_cov))))
sample_cov = sample_cov[!is.infinite(sample_cov)]



mean_len = mean(sample_len)
sd_len = sd(sample_len)

mean_cov = mean(sample_cov)
sd_cov = sd(sample_cov)

sample_len = (sample_len-mean_len)/sd_len

print(sd_cov)
print(mean_cov)


sample_cov = (sample_cov-mean_cov)/sd_cov


data_sample = cbind(sample_len , sample_cov , sampleCon)
data_sample = as.data.frame(data_sample)
colnames(data_sample) = c('sample_len' , 'sample_cov' , 'sampleCon')


print(head(data_sample))


all_len = rep(0 , nrow(all_contacts))
all_cov = rep(0 , nrow(all_contacts))

message('processing all contacts data')
for(i in 1:nrow(all_contacts))
{
  all_len[i] = log(as.numeric(contig_info[contig_info$contig_name==all_contacts[i , 'index1'] , 'length']) * 
                     as.numeric(contig_info[contig_info$contig_name==all_contacts[i , 'index2'] , 'length']))
  
  all_cov[i] = log(as.numeric(contig_info[contig_info$contig_name==all_contacts[i , 'index1'] , 'coverage']) * 
                     as.numeric(contig_info[contig_info$contig_name==all_contacts[i , 'index2'] , 'coverage']))
}

allCon = as.numeric(all_contacts[ , 'contacts'])
all_len = (all_len-mean_len)/sd_len
all_cov = (all_cov-mean_cov)/sd_cov


#print(head(data_sample))

tryCatch(
  {
    message(paste("normalizing",sep=" "))
    
    if(model=='unl'){
      fit1 = glmmTMB(sampleCon~sample_len+sample_cov, data = data_sample,
                     ziformula=~1,family=nbinom2)
      
      
    }else if(model=='nb'){
      fit1 = glmmTMB(sampleCon~sample_len+sample_cov, data = data_sample,
                     ziformula=~0,family=nbinom2)
      
    }else if(model=='none'){
      fit1 = glmmTMB(sampleCon~sample_len+sample_cov, data = data_sample,
                     ziformula=~sample_len+sample_cov , family=nbinom2)
    }
  },
  error = function(e){
    message(e)
    message(paste("\nskip",  sep=" "))
  },
  warning = function(w){
    message(w)
    message(paste("\nskip",  sep=" "))
  }
)

#print("good1")
coeff = as.numeric(fit1$fit$par)
res_sample = sampleCon/exp(coeff[1] + coeff[2]*sample_len + coeff[3]*sample_cov)
mu_sample = exp(coeff[1] + coeff[2]*sample_len + coeff[3]*sample_cov)
index_nonzero = (res_sample > 0)
res_sample_nonzero = res_sample[index_nonzero]
mu_sample_nonzero = mu_sample[index_nonzero]
#print("good2")

res_all =  allCon/exp(coeff[1] + coeff[2]*all_len + coeff[3]*all_cov)
mu_all = exp(coeff[1] + coeff[2]*all_len + coeff[3]*all_cov)
#print(res_all)

#print("good3")
###########detect spurious contacts#################
sigma = summary(fit1)$sigma

pvalue_all = pnbinom(allCon , size = sigma  ,mu = mu_all)
pvalue_sample = pnbinom(sampleCon[index_nonzero] , size =  sigma  , mu = mu_sample_nonzero)

index_spur = (pvalue_all<quantile(pvalue_sample , thres) | res_all < quantile(res_sample_nonzero, thres))


print(index_spur)
print(length(index_spur))

#print(head(all_contacts))


all_contacts$contacts = res_all
#print(all_contacts$contacts)
#print("good4")

res_all_valid = all_contacts[!index_spur, ]
res_all_spur = all_contacts[index_spur, ]

write.table(res_all_valid , file=output_file_name, row.names=F , col.names = F , sep = ',')
#write.table(res_all_spur , file=discarded_contact_name, row.names=F , col.names = F , sep = ',')



