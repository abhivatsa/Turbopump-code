function [red_A,red_B,trans_mat,inv_trans_mat]=serep_modified(A,B)

%     [red_A,red_B,eig_val_A2,eig_vec_right_check]=serep_modified(A,B)
% [eig_val_A2,eig_vec_right_check]=serep_modified(A,B)
    [eig_vec_right,eig_val_A]=eig(A);
    [eig_vec_left,eig_val_A_tran]=eig(A.');
    eig_val_A1=diag(eig_val_A);
    eig_val_Atran1=diag(eig_val_A_tran);
    
    %     Sorting of eigen values and vectors
    [eig_val_A2,eig_val_index]=sort(eig_val_A1);
    [eig_val_Atran2,eig_val_index_tran]=sort(eig_val_Atran1);
    
    for j=1:length(eig_val_A2)
        eig_vec_right1(:,j)=eig_vec_right(:,eig_val_index(j));
        eig_vec_left1(:,j)=eig_vec_left(:,eig_val_index_tran(j));
    end
    
    eig_vec_right_check=eig_vec_right1;
    
    nor_mat=eig_vec_left1.'*eig_vec_right1;
    eig_vec_right1=eig_vec_right1*inv(nor_mat);
    n=length(eig_val_A);

%     node=[1 5 7];
    node=[13 19 26 29 40 43 53];
    
    m_dof=[];
    for cnt=1:length(node)
       m_dof=[m_dof 4*(node(cnt)-1)+1 4*(node(cnt)-1)+2 n/2+4*(node(cnt)-1)+1 n/2+4*(node(cnt)-1)+2];
%        m_dof=[m_dof 4*(node(cnt)-1)+1 4*(node(cnt)-1)+2 4*(node(cnt)-1)+3 4*(node(cnt)-1)+4 n/2+4*(node(cnt)-1)+1 n/2+4*(node(cnt)-1)+2 n/2+4*(node(cnt)-1)+3 n/2+4*(node(cnt)-1)+4];
    end
    m_dof=sort(m_dof);

    mode=1:16;
    
    red1_eig_vec_left=eig_vec_left1(:,mode);
    red1_eig_vec_right=eig_vec_right1(:,mode);
    
    red2_eig_vec_left=red1_eig_vec_left(m_dof,:);
    red2_eig_vec_right=red1_eig_vec_right(m_dof,:);
    
    trans_mat=red2_eig_vec_right*red1_eig_vec_left.';
    inv_trans_mat=red1_eig_vec_right*red2_eig_vec_left.'*pinv(red2_eig_vec_left.')*pinv(red2_eig_vec_right);
    
    red_A=trans_mat*A*inv_trans_mat;
    red_B=trans_mat*B;
%     check=sort(eig(red_A));
end
