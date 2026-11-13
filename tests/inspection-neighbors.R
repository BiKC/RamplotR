source('shinyRam/R/inspection.R')
set.seed(28)
target <- data.frame(chain='A',resno=1L,resid='ALA',insert='',
 elety=c('N','CA','C'),elesy=c('N','C','C'),x=c(0,1,2),y=0,z=0)
hetero <- data.frame(chain='Z',resno=rep(1:100,each=2),resid='LIG',insert='',
 elety=rep(c('C1','C2'),100),elesy='C',x=runif(200,-10,10),
 y=runif(200,-10,10),z=runif(200,-10,10))
hetero[1:2,c('x','y','z')] <- 0
pdb <- list(atom=target,hetero_atom=hetero)
row <- data.frame(chain='A',resi=1L,resn='ALA',insertion_code='')
# Independent all-pairs reference; production uses bounding-box rejection.
distance <- sqrt(outer(hetero$x,target$x,'-')^2+
 outer(hetero$y,target$y,'-')^2+outer(hetero$z,target$z,'-')^2)
nearest <- apply(distance,1,min)
for (radius in c(0,2,6,50)) {
 got <- ram_nearby_hetero_context(pdb,row,max_distance=radius,max_hits=100L)
 keep <- which(nearest<=radius)
 best <- vapply(split(keep,hetero$resno[keep]),function(ix)ix[which.min(nearest[ix])],integer(1))
 best <- best[order(nearest[best],hetero$resno[best])]
 stopifnot(identical(as.integer(got$resi),as.integer(hetero$resno[best])),
  isTRUE(all.equal(got$distance,unname(nearest[best]),tolerance=1e-12)),
  identical(got$target_atom,target$elety[apply(distance[best,,drop=FALSE],1,which.min)]))
}
stopifnot(nrow(ram_nearby_hetero_context(pdb,row,max_distance=-1))==0L)
message('Vectorized local context matches independent all-pairs distances and nearest atom identities.')
