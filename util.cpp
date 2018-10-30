//
// Created by andyshen on 18-10-25.
//

#include<iostream>
#include <vector>
#include <algorithm>
#define DEBUG 0
float calutate_acc(std::vector<int>lsh, std::vector<int> v_lsh){
    sort(lsh.begin(), lsh.end());
    sort(v_lsh.begin(), v_lsh.end());
    float count = 0.0;
    for(std::vector<int>::const_iterator lsh_index=lsh.cbegin(); lsh_index!=lsh.cend(); lsh_index++){
        for(std::vector<int>::const_iterator vlsh_index=v_lsh.cbegin(); vlsh_index!=v_lsh.cend(); vlsh_index++){
            if(*lsh_index<*vlsh_index){
                break;
            }
            else if(*lsh_index==*vlsh_index){
                count++;
            }
        }
    }
    std::cout<<"sim: "<<count<<std::endl;
    return count/lsh.size();
}
#if DEBUG
int main(){
    std::vector<int> a = {1,3,8,4,5,6,7};
    std::vector<int> b = {1,2,8,10,5,6,7};
    std::cout<<calutate_acc(a,b)<<std::endl;
    return 0;
}
#endif
