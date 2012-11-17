%%%
%%% Copyright 2011, Boundary
%%%
%%% Licensed under the Apache License, Version 2.0 (the "License");
%%% you may not use this file except in compliance with the License.
%%% You may obtain a copy of the License at
%%%
%%%     http://www.apache.org/licenses/LICENSE-2.0
%%%
%%% Unless required by applicable law or agreed to in writing, software
%%% distributed under the License is distributed on an "AS IS" BASIS,
%%% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%%% See the License for the specific language governing permissions and
%%% limitations under the License.
%%%


%%%-------------------------------------------------------------------
%%% File:      bear.erl
%%% @author    joe williams <j@boundary.com>
%%% @doc
%%% statistics functions for calucating based on id and a list of values
%%% @end
%%%------------------------------------------------------------------

-module(bear).

-compile([export_all]).

-export([get_statistics/1, get_statistics/2, get_correlation_statistics/2]).

-record(scan_result, {
    n=0,
    sumX=0,
    sumXX=0,
    sumInv=0,
    sumLog,
    max,
    min,
    x2,
    x3,
    x4
}).

-compile([native]).

get_statistics(Values) ->
    DefaultStats = [
        min,
        max,
        n,
        arithmetic_mean,
        geometric_mean,
        harmonic_mean,
        variance,
        standard_deviation,
        skewness,
        kurtosis,
        {percentile, [0.5, 0.75, 0.99, 0.999]},
        histogram
    ],
    get_statistics(Values, DefaultStats).

get_statistics(Values, Stats) ->
    ScanRes = scan_values(Values),
    compute_statistics(Values, Stats, ScanRes, []).

get_correlation_statistics(V1, V2) when length(V1) == length(V2) ->
    [
     {covariance, get_covariance(V1, V2)},
     {tau, get_kendall_correlation(V1, V2)},
     {rho, get_pearson_correlation(V1, V2)},
     {r, get_spearman_correlation(V1, V2)}
    ].

%%%===================================================================
%%% Internal functions
%%%===================================================================

compute_statistics(Values, [min|S], SR, Acc) ->
    compute_statistics(Values, S, SR, [{min, SR#scan_result.min}|Acc]);
compute_statistics(Values, [max|S], SR, Acc) ->
    compute_statistics(Values, S, SR, [{max, SR#scan_result.max}|Acc]);
compute_statistics(Values, [n|S], SR, Acc) ->
    compute_statistics(Values, S, SR, [{n, SR#scan_result.n}|Acc]);
compute_statistics(Values, [arithmetic_mean|S], SR, Acc) ->
    Stat = arithmetic_mean(SR),
    compute_statistics(Values, S, SR, [{arithmetic_mean, Stat}|Acc]);
compute_statistics(Values, [geometric_mean|S], SR, Acc) ->
    Stat = geometric_mean(SR),
    compute_statistics(Values, S, SR, [{geometric_mean, Stat}|Acc]);
compute_statistics(Values, [harmonic_mean|S], SR, Acc) ->
    Stat = harmonic_mean(SR),
    compute_statistics(Values, S, SR, [{harmonic_mean, Stat}|Acc]);
compute_statistics(Values, [{percentile, Stats}|S], SR, Acc) ->
    Sorted = lists:sort(Values),
    Percentiles = [{Stat, percentile(Sorted, SR, Stat)} || Stat <- Stats],
    compute_statistics(Values, S, SR, [{percentile, Percentiles}|Acc]);
compute_statistics(Values, [variance|S], #scan_result{x2=undefined}=SR, Acc) ->
    compute_statistics(Values, [variance|S], scan_values2(Values, SR), Acc);
compute_statistics(Values, [variance|S], SR, Acc) ->
    Stat = variance(Values, SR),
    compute_statistics(Values, S, SR, [{variance, Stat}|Acc]);
compute_statistics(Values, [standard_deviation|S], #scan_result{x2=undefined}=SR, Acc) ->
    compute_statistics(Values, [standard_deviation|S], scan_values2(Values, SR), Acc);
compute_statistics(Values, [standard_deviation|S], SR, Acc) ->
    Stat = standard_deviation(Values, SR),
    compute_statistics(Values, S, SR, [{standard_deviation, Stat}|Acc]);
compute_statistics(Values, [skewness|S], #scan_result{x3=undefined}=SR, Acc) ->
    compute_statistics(Values, [skewness|S], scan_values2(Values, SR), Acc);
compute_statistics(Values, [skewness|S], SR, Acc) ->
    Stat = skewness(Values, SR),
    compute_statistics(Values, S, SR, [{skewness, Stat}|Acc]);
compute_statistics(Values, [kurtosis|S], #scan_result{x4=undefined}=SR, Acc) ->
    compute_statistics(Values, [kurtosis|S], scan_values2(Values, SR), Acc);
compute_statistics(Values, [kurtosis|S], SR, Acc) ->
    Stat = kurtosis(Values, SR),
    compute_statistics(Values, S, SR, [{kurtosis, Stat}|Acc]);
compute_statistics(Values, [histogram|S], #scan_result{x2=undefined}=SR, Acc) ->
    compute_statistics(Values, [histogram|S], scan_values2(Values, SR), Acc);
compute_statistics(Values, [histogram|S], SR, Acc) ->
    Stat = histogram(Values, SR),
    compute_statistics(Values, S, SR, [{histogram, Stat}|Acc]);
compute_statistics(_, [_|_], _, _) ->
    {error, unknown_statistic};
compute_statistics(_, [], _, Acc) ->
    Acc.

scan_values([X|Values]) ->
    scan_values(Values, #scan_result{n=1, sumX=X, sumXX=X*X,
             sumLog=math_log(X),
             max=X, min=X, sumInv=inverse(X)}).

scan_values([X|Values],
      #scan_result{n=N, sumX=SumX, sumXX=SumXX, sumLog=SumLog,
                   max=Max, min=Min, sumInv=SumInv}=Acc) ->
    scan_values(Values,
                Acc#scan_result{n=N+1, sumX=SumX+X, sumXX=SumXX+X*X,
                                sumLog=SumLog+math_log(X),
                                max=max(X,Max), min=min(X,Min),
                                sumInv=SumInv+inverse(X)});
scan_values([], Acc) ->
    Acc.

scan_values2(Values, #scan_result{n=N, sumX=SumX}=SR0) ->
    SR1 = SR0#scan_result{x2=0, x3=0, x4=0},
    scan_values2(Values, SumX/N, SR1).

scan_values2([X|Values], Mean, #scan_result{x2=X2, x3=X3, x4=X4}=Acc) ->
    Diff = X-Mean,
    Diff2 = Diff*Diff,
    Diff3 = Diff2*Diff,
    Diff4 = Diff2*Diff2,
    scan_values2(Values, Mean, Acc#scan_result{x2=X2+Diff2, x3=X3+Diff3,
            x4=X4+Diff4});
scan_values2([], _, Acc) ->
    Acc.

arithmetic_mean(#scan_result{n=N, sumX=Sum}) ->
    Sum/N.

geometric_mean(#scan_result{n=N, sumLog=SumLog}) ->
    math:exp(SumLog/N).

harmonic_mean(#scan_result{n=N, sumInv=Sum}) ->
    N/Sum.

percentile(SortedValues, #scan_result{n=N}, Percentile)
  when is_list(SortedValues) ->
    Element = round(Percentile * N),
    lists:nth(Element, SortedValues).

%% Two pass variance
%% Results match those given by the 'var' function in R
variance(_, #scan_result{n=N, x2=X2}) ->
    X2/(N-1).

standard_deviation(Values, ScanRes) ->
    math:sqrt(variance(Values, ScanRes)).

%% http://en.wikipedia.org/wiki/Skewness
%%
%% skewness results should match this R function:
%% skewness <- function(x) {
%%    m3 <- mean((x - mean(x))^3)
%%    skew <- m3 / (sd(x)^3)
%%    skew
%% }
skewness(Values, #scan_result{n=N, x3=X3}=ScanRes) ->
    case math:pow(standard_deviation(Values, ScanRes), 3) of
        0.0 ->
            0.0;  %% Is this really the correct thing to do here?
        Else ->
            (X3/N)/Else
    end.

%% http://en.wikipedia.org/wiki/Kurtosis
%%
%% results should match this R function:
%% kurtosis <- function(x) {
%%     m4 <- mean((x - mean(x))^4)
%%     kurt <- m4 / (sd(x)^4) - 3
%%     kurt
%% }
kurtosis(Values, #scan_result{n=N, x4=X4}=ScanRes) ->
    case math:pow(standard_deviation(Values, ScanRes), 4) of
        0.0 ->
            0.0;  %% Is this really the correct thing to do here?
        Else ->
            ((X4/N)/Else) - 3
    end.

histogram(Values, ScanRes) ->
    Bins = get_hist_bins(ScanRes#scan_result.min,
                         ScanRes#scan_result.max,
                         standard_deviation(Values, ScanRes),
                         length(Values)
                        ),

    Dict = lists:foldl(fun (Value, Dict) ->
             update_bin(Value, Bins, Dict)
           end,
           dict:from_list([{Bin, 0} || Bin <- Bins]),
           Values),

    lists:sort(dict:to_list(Dict)).

update_bin(Value, [Bin|_Bins], Dict) when Value =< Bin ->
    dict:update_counter(Bin, 1, Dict);
update_bin(Values, [_Bin|Bins], Dict) ->
    update_bin(Values, Bins, Dict).

%% two pass covariance
%% (http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Covariance)
%% matches results given by excel's 'covar' function
get_covariance(Values1, Values2) ->
    {SumX, SumY, N} = foldl2(fun (X, Y, {SumX, SumY, N}) ->
             {SumX+X, SumY+Y, N+1}
           end, {0,0,0}, Values1, Values2),
    MeanX = SumX/N,
    MeanY = SumY/N,
    Sum = foldl2(fun (X, Y, Sum) ->
       Sum + ((X - MeanX) * (Y - MeanY))
     end,
     0, Values1, Values2),
    Sum/N.

get_kendall_correlation(Values1, Values2) ->
    bear:kendall_correlation(Values1, Values2).

get_spearman_correlation(Values1, Values2) ->
    TR1 = ranks_of(Values1),
    TR2 = ranks_of(Values2),
    Numerator   = 6 * foldl2(fun (X, Y, Acc) ->
             Diff = X-Y,
             Acc + Diff*Diff
           end, 0, TR1,TR2),
    N = length(Values1),
    Denominator = math:pow(N,3)-N,
    1-(Numerator/Denominator).

ranks_of(Values) when is_list(Values) ->
    [Fst|Rest] = revsort(Values),
    TRs = ranks_of(Rest, [], 2, Fst, 1),
    Dict = gb_trees:from_orddict(TRs),
    L = lists:foldl(fun (Val, Acc) ->
                            Rank = gb_trees:get(Val, Dict),
                            [Rank|Acc]
                    end, [], Values),
    lists:reverse(L).

ranks_of([E|Es],Acc, N, E, S) ->
    ranks_of(Es, Acc, N+1, E, S);
ranks_of([E|Es], Acc, N, P, S) ->
    ranks_of(Es,[{P,(S+N-1)/2}|Acc], N+1, E, N);
ranks_of([],  Acc, N, P, S) ->
    [{P,(S+N-1)/2}|Acc].

get_pearson_correlation(Values1, Values2) ->
    {SumX, SumY, SumXX, SumYY, SumXY, N} =
  foldl2(fun (X,Y,{SX, SY, SXX, SYY, SXY, N}) ->
          {SX+X, SY+Y, SXX+X*X, SYY+Y*Y, SXY+X*Y, N+1}
        end, {0,0,0,0,0,0}, Values1, Values2),
    Numer = (N*SumXY) - (SumX * SumY),
    case math:sqrt(((N*SumXX)-(SumX*SumX)) * ((N*SumYY)-(SumY*SumY))) of
        0.0 ->
            0.0; %% Is this really the correct thing to do here?
        Denom ->
            Numer/Denom
    end.

revsort(L) ->
    lists:reverse(lists:sort(L)).

%% Foldl over two lists
foldl2(F, Acc, [I1|L1], [I2|L2]) when is_function(F,3) ->
    foldl2(F, F(I1, I2, Acc), L1, L2);
foldl2(_F, Acc, [], []) ->
    Acc.

%% wrapper for math:log/1 to avoid dividing by zero
math_log(0) ->
    1;
math_log(X) ->
    math:log(X).

%% wrapper for calculating inverse to avoid dividing by zero
inverse(0) ->
    0;
inverse(X) ->
    1/X.

get_hist_bins(Min, Max, StdDev, Count) ->
    BinWidth = get_bin_width(StdDev, Count),
    BinCount = get_bin_count(Min, Max, BinWidth),
    case get_bin_list(BinWidth, BinCount, []) of
        List when length(List) =< 1 ->
            [Max];
        Bins ->
            %% add Min to Bins
            [Bin + Min || Bin <- Bins]
    end.

get_bin_list(Width, Bins, Acc) when Bins > length(Acc) ->
    Bin = ((length(Acc) + 1) * Width ),
    get_bin_list(Width, Bins, [round_bin(Bin)| Acc]);
get_bin_list(_, _, Acc) ->
    lists:usort(Acc).

round_bin(Bin) ->
    Base = case erlang:trunc(math:pow(10, round(math:log10(Bin) - 1))) of
        0 ->
            1;
        Else ->
            Else
    end,
    %io:format("bin ~p, base ~p~n", [Bin, Base]),
    round_bin(Bin, Base).

round_bin(Bin, Base) when Bin rem Base == 0 ->
    Bin;
round_bin(Bin, Base) ->
    Bin + Base - (Bin rem Base).

% the following is up for debate as far as what the best method
% of choosing bin counts and widths. these seem to work *good enough*
% in my testing

% bin width based on Sturges
% http://www.jstor.org/pss/2965501
get_bin_width(StdDev, Count) ->
    %io:format("stddev: ~p, count: ~p~n", [StdDev, Count]),
    case round((3.5 * StdDev) / math:pow(Count, 0.3333333)) of
        0 ->
            1;
        Else ->
            Else
    end.

% based on the simple ceilng function at
% http://en.wikipedia.org/wiki/Histograms#Number_of_bins_and_width
% with a modification to attempt to get on bin beyond the max value
get_bin_count(Min, Max, Width) ->
    %io:format("min: ~p, max: ~p, width ~p~n", [Min, Max, Width]),
    round((Max - Min) / Width) + 1.

%% taken from http://crunchyd.com/scutil/
%% All code here is MIT Licensed
%%  http://scutil.com/license.html

% seems to match the value returned by the 'cor' (method="kendal") R function
% http://en.wikipedia.org/wiki/Kendall_tau_rank_correlation_coefficient
kendall_correlation(List1, List2) when is_list(List1), is_list(List2) ->
    {RA,_} = lists:unzip(tied_ordered_ranking(List1)),
    {RB,_} = lists:unzip(tied_ordered_ranking(List2)),

    Ordering = lists:keysort(1, lists:zip(RA,RB)),
    {_,OrdB} = lists:unzip(Ordering),

    N = length(List1),
    P = lists:sum(kendall_right_of(OrdB, [])),

    -(( (4*P) / (N * (N - 1))) - 1).

simple_ranking(List) when is_list(List) ->
    lists:zip(lists:seq(1,length(List)),lists:reverse(lists:sort(List))).

tied_ranking(List) ->
    tied_rank_worker(simple_ranking(List), [], no_prev_value).

tied_ordered_ranking(List) when is_list(List) ->
    tied_ordered_ranking(List, tied_ranking(List), []).

tied_ordered_ranking([], [], Work) ->
    lists:reverse(Work);

tied_ordered_ranking([Front|Rem], Ranks, Work) ->
    {value,Item}  = lists:keysearch(Front,2,Ranks),
    {IRank,Front} = Item,
    tied_ordered_ranking(Rem, Ranks--[Item], [{IRank,Front}]++Work).

kendall_right_of([], Work) ->
    lists:reverse(Work);
kendall_right_of([F|R], Work) ->
    kendall_right_of(R, [kendall_right_of_item(F,R)]++Work).

kendall_right_of_item(B, Rem) ->
    length([R || R <- Rem, R < B]).

tied_add_prev(Work, {FoundAt, NewValue}) ->
    lists:duplicate( length(FoundAt), {lists:sum(FoundAt)/length(FoundAt), NewValue} ) ++ Work.

tied_rank_worker([], Work, PrevValue) ->
    lists:reverse(tied_add_prev(Work, PrevValue));

tied_rank_worker([Item|Remainder], Work, PrevValue) ->
    case PrevValue of
        no_prev_value ->
            {BaseRank,BaseVal} = Item,
            tied_rank_worker(Remainder, Work, {[BaseRank],BaseVal});
        {FoundAt,OldVal} ->
            case Item of
                {Id,OldVal} ->
                    tied_rank_worker(Remainder, Work, {[Id]++FoundAt,OldVal});
                {Id,NewVal} ->
                    tied_rank_worker(Remainder, tied_add_prev(Work, PrevValue), {[Id],NewVal})

            end
    end.

