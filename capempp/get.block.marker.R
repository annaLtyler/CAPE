get.block.marker <- function(data.obj, block.name, collapsed.net){

	if(collapsed.net){
		return(data.obj$linkage.blocks.collapsed[which(names(data.obj$linkage.blocks.collapsed) == block.name)])
		}else{
		return(data.obj$linkage.blocks.full[which(names(data.obj$linkage.blocks.full) == block.name)])			
		}
}
