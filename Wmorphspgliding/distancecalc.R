distancecalc <- function(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10){
  distance <- data.frame(matrix(nrow = nrow(pc_results_labman), ncol = 3))
  names(distance) <- c("dis_all10","dis_first6","dis_first2")
  
  output <- data.frame(matrix(nrow = 1, ncol = 6))
  names(output) <- c("closest_first2","closest_first6","closest_all10","mindis_first2","mindis_first6","mindis_all10")
  
  for (i in 1:nrow(pc_results_labman)){
  d_pc1 = (pc_results_labman$PC1[i]-PC1)^2
  d_pc2 = (pc_results_labman$PC2[i]-PC2)^2
  d_pc3 = (pc_results_labman$PC3[i]-PC3)^2
  d_pc4 = (pc_results_labman$PC4[i]-PC4)^2
  d_pc5 = (pc_results_labman$PC5[i]-PC5)^2
  d_pc6 = (pc_results_labman$PC6[i]-PC6)^2
  d_pc7 = (pc_results_labman$PC7[i]-PC7)^2
  d_pc8 = (pc_results_labman$PC8[i]-PC8)^2
  d_pc9 = (pc_results_labman$PC9[i]-PC9)^2
  d_pc10 = (pc_results_labman$PC10[i]-PC10)^2
  
  distance$dis_all10[i]  = sqrt(d_pc1+d_pc2+d_pc3+d_pc4+d_pc5+d_pc6+d_pc7+d_pc8+d_pc9+d_pc10)
  distance$dis_first6[i] = sqrt(d_pc1+d_pc2+d_pc3+d_pc4+d_pc5+d_pc6)
  distance$dis_first2[i] = sqrt(d_pc1+d_pc2)
  }
  
  output$closest_first2 = pc_results_labman$elbow.angle[which.min(distance$dis_first2)]
  output$closest_first6 = pc_results_labman$elbow.angle[which.min(distance$dis_first6)]
  output$closest_all10  = pc_results_labman$elbow.angle[which.min(distance$dis_all10)]
  
  output$mindis_first2 = distance$dis_first2[which.min(distance$dis_first2)]
  output$mindis_first6 = distance$dis_first6[which.min(distance$dis_first6)]
  output$mindis_all10  = distance$dis_all10[which.min(distance$dis_all10)]
  
  tmp <- pc_results_labman[order(distance$dis_first2),]
  tmp2 <- distance[order(distance$dis_first2),]
  output$closest_first2_avg3 = (tmp$elbow.angle[1]*tmp2$dis_first2[1]+tmp$elbow.angle[2]*tmp2$dis_first2[2]+tmp$elbow.angle[3]*tmp2$dis_first2[3])/(tmp2$dis_first2[1]+tmp2$dis_first2[2]+tmp2$dis_first2[3])
  output$mindis_first2_avg3  = (tmp2$dis_first2[1]+tmp2$dis_first2[2]+tmp2$dis_first2[3])/3
  
  tmp <- pc_results_labman[order(distance$dis_first6),]
  tmp2 <- distance[order(distance$dis_first6),]
  output$closest_first6_avg3 = (tmp$elbow.angle[1]*tmp2$dis_first2[1]+tmp$elbow.angle[2]*tmp2$dis_first2[2]+tmp$elbow.angle[3]*tmp2$dis_first2[3])/(tmp2$dis_first2[1]+tmp2$dis_first2[2]+tmp2$dis_first2[3])
  output$mindis_first6_avg3  = (tmp2$dis_first2[1]+tmp2$dis_first2[2]+tmp2$dis_first2[3])/3
  
  tmp <- pc_results_labman[order(distance$dis_all10),]
  tmp2 <- distance[order(distance$dis_all10),]
  output$closest_all10_avg3 = (tmp$elbow.angle[1]*tmp2$dis_first2[1]+tmp$elbow.angle[2]*tmp2$dis_first2[2]+tmp$elbow.angle[3]*tmp2$dis_first2[3])/(tmp2$dis_first2[1]+tmp2$dis_first2[2]+tmp2$dis_first2[3])
  output$mindis_all10_avg3  = (tmp2$dis_first2[1]+tmp2$dis_first2[2]+tmp2$dis_first2[3])/3
  
  return(output)
}