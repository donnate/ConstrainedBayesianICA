library(igraph)
library(rmutil)


make_small_subnetwork <- function(n0, GG="NA", eta = 0.5){
    if (GG =="NA"){
        GG <- sample_pa(n0, power=2, m=1, directed=FALSE)
    }

  #Sigma_i <-  matrix(0,N,N)
  #Sigma_i[a:(a+14), a:(a+14)] <- as_adjacency_matrix(GG, sparse = FALSE)
  Sigma_i = as_adjacency_matrix(GG, sparse = FALSE) 
  edges <- runif(sum(Sigma_i)/2, 0.2, 0.4) 
  Sigma_i[lower.tri(Sigma_i,diag = FALSE)][Sigma_i[lower.tri(Sigma_i,diag = FALSE)]>0] = edges
  Sigma_i = Sigma_i + t(Sigma_i) #### weighted adjacency matrix
  L = diag(apply(Sigma_i, 1,sum)) - Sigma_i + eta * diag(n0)  ### Rgularized Laplacian for the graph structure
  svd_Om = svd(L)
  Sigma = svd_Om$u %*% diag(sapply(svd_Om$d, FUN = function(x){ifelse(x>10^(-7), 1/x, 0)})) %*% t(svd_Om$u)
  diag_Sigma = diag(sqrt(1.0/diag(Sigma)))
  Sigma = 1/sqrt(n0) * diag_Sigma  %*% Sigma  %*% diag_Sigma 
  return(list(Sigma = Sigma, Adj=Sigma_i))
}


create_graph <- function(N,K, n0){
  Sigma_i <- list()
  Adj = matrix(0, N,N)
  a = 1
  list_links <- c()
  for (i in 1:K){
    test  = make_small_subnetwork(n0[i], GG="NA")
    Sigma_i[[i]] <-  test$Sigma
    if (i==1){
        indices = a:(a+n0[i]-1)
    }
    else{
        indices = c(list_links, a:(a+n0[i]-length(list_links )-1))
    }
    #print("list indices")
    print(indices)
    #print(length(indices))
    Adj[indices, indices] = test$Adj
    a = a+n0[i]-length(list_links)
    #print(a)
    #print("list links")
    list_links <- c(list_links,a-1)
    print(list_links)
  }
  Sigma_i[[K+1]] <- diag(N - a + 1)
  return(list(Sigma_i= Sigma_i, Adj = Adj, Graph =graph_from_adjacency_matrix(Adj)))
}


create_connected_graph <- function(N,K, n0){
  Sigma_i <- list()
  Adj = matrix(0, N,N)
  a = 1
  list_links <- c()
  for (i in 1:K){
    test  = make_small_subnetwork(n0[i], GG="NA")
    Sigma_i[[i]] <-  test$Sigma
    if (i==1){
      indices = a:(a+n0[i]-1)
    }
    else{
      indices = c(list_links, a:(a+n0[i]-length(list_links )-1))
    }
    #print("list indices")

    #print(length(indices))
    Adj[indices, indices] = test$Adj
    a = a+n0[i]-length(list_links)
    #print(a)
    #print("list links")
    list_links <- c(list_links,a-1)
    print(list_links)
  }
  Adj = Adj[1:(a-1), 1:(a-1)]
  Graph = graph_from_adjacency_matrix(Adj)
  for (v in 1:length(V(Graph))){
    vertex_attr(Graph, "name", v) = v
  }
  return(list(Sigma_i= Sigma_i, Adj = Adj, Graph =Graph))
}



create_block_factor_model <- function(N,K, n0, Time, mu){
  graphs <-  create_graph(N,K, n0)
  A = matrix(0, K, N)
  a = 1
  for (i in 1:K){
    A[i, a:(a+n0[i]-1)] = mvrnorm(1, rep(mu[i], n0[i]), graphs$Sigma_i[[i]])
    a = a + n0[i]
  }
  
  return(A)
}

graph2cov <- function(Adj, eta = 0.5){
  N = nrow(Adj)
  L_structural = diag(apply(Adj, 1,sum)) -Adj + eta* diag(N)  ### Rgularized Laplacian for the graph structure
  svd_struct = svd(L_structural)
  Sigma = svd_struct$u %*% diag(1/svd_struct$d) %*% t(svd_struct$u)
  D_S = diag(sqrt(1/diag(Sigma)))
  Sigma = D_S %*% Sigma %*% D_S
  Omega=svd(Sigma)
  return(list(Sigma=Sigma, Omega=Omega$u %*% diag(1/Omega$d) %*% t(Omega$u)))
}


scramble <- function(Adj, Graph, noise_level, n0, eta = 0.5){
  K = length(n0)
  N = nrow(Adj)
  n_edges = sum(Adj!=0)/2
  n_edges2change = ceiling(n_edges *  noise_level)
  list_e = get.edgelist(Graph)
  new_Adj = Adj
  list_e= list_e[list_e[,1]<list_e[,2],]
  list_edges2change = sample(1:n_edges,n_edges2change, replace = FALSE )
  for (i in list_edges2change){
    a = list_e[i,1]
    b = list_e[i,2]
    candidates_b = setdiff(neighbors(Graph,b), neighbors(Graph,a))
    if (length(candidates_b) == 1 && candidates_b[1] == a){
      ### choose neighbor c from a and connect b and c
      candidates_a = setdiff(neighbors(Graph,a), c(b,neighbors(Graph,b)))
      c = sample(candidates_a,1)
      new_Adj[b,c] = runif(1,0.2,0.4)
      new_Adj[c,b] = new_Adj[b,c]
    }
    else{
      #### choose another neighbor of b
      c = sample(candidates_b,1)
      new_Adj[a,c] = Adj[a,b]
      new_Adj[c,a] = new_Adj[a,c]
      new_Adj[a,b] = 0
      new_Adj[b,a] = 0 
    }

    
  }
  
  
  G=graph_from_adjacency_matrix(new_Adj)
  L_structural = diag(apply(new_Adj, 1,sum)) - new_Adj + eta* diag(N)  ### Rgularized Laplacian for the graph structure
  svd_struct = svd(L_structural)
  Sigma = svd_struct$u %*% diag(1/svd_struct$d) %*% t(svd_struct$u)
  D_S = diag(sqrt(1/diag(Sigma)))
  Sigma = D_S %*% Sigma %*% D_S
  Omega=svd(Sigma)
  
  
  Sigma_i <- list()
  a = 1
  list_links <- c()
  for (i in 1:K){
      if (i==1){
          indices = a:(a+n0[i]-1)
      }
      else{
          indices = c(list_links, a:(a+n0[i]-length(list_links )-1))
      }
      print(indices)
      Sigma_i[[i]] <- Sigma[indices, indices]
      a = a+n0[i]-length(list_links)
      list_links <- c(list_links,a-1)
  }
  Sigma_i[[K+1]] <- Sigma[a:N, a:N]
  return(list(Adj =as_adjacency_matrix(G, sparse = FALSE),
              G=G, Omega=Omega, Sigma_i = Sigma_i))
}


scramble_connected <- function(Graph, nodes2swap){
  GG = graph_from_adjacency_matrix(Graph$Adj)
  for (v in 1:length(V(GG))){
    vertex_attr(GG, "name", v) = v
  }
  for(i in 1:(length(nodes2swap)/2)){
      vertex_attr(GG,"name",nodes2swap[2*(i-1)+ 1] ) = nodes2swap[2*i]
      vertex_attr(GG,"name",nodes2swap[2*i] ) = nodes2swap[2*(i-1)+ 1]
    }
    return(list(Adj =Graph$Adj,
    Graph=GG, Omega=Graph$Omega, Sigma_i = Graph$Sigma_i))
}


