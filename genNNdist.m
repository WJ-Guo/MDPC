function [pd,affinity,distNN,affinityNN,kneigor] = genNNdist(data,k_nNum,sigma)
%%dist,affinity-----distance matrix,and affinity matrix
%%distNN,affinityNN------k_nNum-nearest neighborhood distance matrix and affinity
% matrix,
% data ---NxM datasets

			
			NT = size(data,1);
			distNN= ones(NT,NT)*inf;
			affinityNN=zeros(NT,NT);
			kneigor = ones(NT,k_nNum)*inf;
			%这是个非常重要的变量，它可以使得距离矩阵变得非常稀疏


			pd = dist(data');
			sigma = mean(mean(pd));
			
			affinity = exp(-(pd.^2)/(2*sigma^2));
			for i=1:NT    
				%for j=1:NT
				%    dist(i,j) = sqrt((data(i,1) - data(j,1))^2 + (data(i,2) - data(j,2))^2); 
				%    affinity(i,j) = exp(-dist(i,j)^2/(2*sigma^2));
				%end
								d=pd(i,:);
								d(i)=Inf;%因为d中自身距离是零，这个不能算
						if(k_nNum<log(NT))
						   min_id=[];
							for k=1:k_nNum
							  [tmp,id]=min(d);
							  d(id)=Inf;
							  min_id=[min_id id];
							  % sort>=O(n*logn),so we take min: O(n).total time:O(k*n)
							end
						else
						  [tmp,id]=sort(d);
						  min_id=id(1:k_nNum);%最近邻下标
						end
						distNN(i,min_id) = pd(i,min_id);
						kneigor(i,:) = min_id;
						affinityNN(i,min_id) = exp(-pd(i,min_id).^2./(2*sigma^2));
			end
			
            
			disp('Computing symmetric affinity matrix...')
			A1 = triu(affinityNN);
			A1 = A1 + A1';
			A2 = tril(affinityNN);
			A2 = A2 + A2';
			clear affinityNN;
			
			%互近邻集合
			affinityNN = max(A1, A2);
            affinityNN(1:NT+1:end) = 1;
			
			disp('Computing symmetric distance matrix...')
			%distNN(1,find(distNN(1,:)~=Inf))
			%distNN(find(distNN(1,:)~=Inf),1)
			%pause
			
			
			
			A1 = triu(distNN);
			A1 = A1 + A1';
			%A1(1,find(A1(1,:)~=Inf))
			%A1(find(A1(1,:)~=Inf),1)
			%pause

			A2 = tril(distNN);
			A2 = A2 + A2';
			%A2(1,find(A2(1,:)~=Inf))
			%A2(find(A2(1,:)~=Inf),1)
			%pause
			clear distNN;
			distNN = max(A1, A2);
            distNN(1:NT+1:end) = 0;