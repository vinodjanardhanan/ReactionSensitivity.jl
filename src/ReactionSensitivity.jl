module ReactionSensitivity
using LightXML, Printf
using DifferentialEquations, Statistics, GlobalSensitivity, QuasiMonteCarlo, Sundials
using ReactionCommons, RxnHelperUtils, SurfaceReactions
using BatchReactor, PlugFlowReactor, StirredReactor

export  rxn_gsa
global o_streams
global gsa_model

struct Solver 
    abstol::Float64
    reltol::Float64
end

struct  Gsa_parameters 
    monitor::Array{String,1}   
    lower_pc::Float64
    upper_pc::Float64
    N::Int64
    gsa_srxn_ids::Array{Int64,1}
    gsa_srxn_constraint_ids::Array{Int64,1}    
end
struct Mechanisms 
    smd::MechanismDefinition
end


function rxn_gsa(input_file::AbstractString, lib_dir::AbstractString)
    xmldoc = parse_file(input_file)
    xmlroot = root(xmldoc)

    

    # list of gasphase species 
    gsa_species = get_collection_from_xml(xmlroot,"gsa_species")
    # lower bound for parameters (in percentages)
    lower_pc = get_value_from_xml(xmlroot,"gsa_lb")
    # upper bound for parameters (in percentages)
    upper_pc = get_value_from_xml(xmlroot,"gsa_ub")
    
    # Size of the sample
    sample_N = Int(get_value_from_xml(xmlroot,"gsa_N"))

    # get the reactions id which forms the parameter 
    gsa_srxn_ids,gsa_srxn_constraint_ids = get_rxn_ids(xmlroot,"p_smech")              
    gsa_p = Gsa_parameters(gsa_species,lower_pc,upper_pc,sample_N,gsa_srxn_ids,gsa_srxn_constraint_ids)
    

    # get the reactor model 
    model = get_text_from_xml(xmlroot,"gsa_model")

    if Sys.islinux() || Sys.isapple()
        dirs = split(input_file, "/")        
    else
        dirs = split(input_file, "\\")
    end
    pop!(dirs)   
    push!(dirs, model)             
    model = joinpath(dirs)

    
    model_root = root(parse_file(model))
    global gsa_model = name(model_root)


    abstol = get_value_from_xml(xmlroot,"abstol")
    reltol = get_value_from_xml(xmlroot, "reltol")
    soln_cntrl = Solver(abstol, reltol)

    #= output file generation =#
    # header_data = [ "p"*string(k) for k in gsa_p.gsa_srxn_ids]    
    # println(header_data)
    sobol_total = open(output_file(input_file, "sobol_total.dat"),"w")     
    sobol_first = open(output_file(input_file, "sobol_first.dat"),"w")    
    global o_streams = sobol_first, sobol_total
    create_header(sobol_total,["parameter"], gsa_p.monitor)
    create_header(sobol_first,["parameter"], gsa_p.monitor)


    if lowercase(strip(gsa_model)) == "batch"
        println("Performing global sensitivity analysis using Batch reactor model")
        chem = Chemistry(true, false, false, f->())
        params, prob, t_span = batch_reactor(model,lib_dir,true, chem)        
        retcode = gsa_problem(params, prob, t_span, gsa_p, soln_cntrl)       
        return retcode
    elseif lowercase(strip(gsa_model)) == "cstr"
        println("Performing global sensitivity analysis using CSTR model")
        chem = Chemistry(true, false, false, f->())
        params, prob, t_span = cstr(model,lib_dir,true, chem)        
        retcode = gsa_problem(params, prob, t_span, gsa_p, soln_cntrl)               
        return retcode
    elseif lowercase(strip(gsa_model)) == "plug"
        println("Performing global sensitivity analysis using Plug flow model")
        chem = Chemistry(true, false, false, f->())
        params, prob, t_span = plug(model,lib_dir,true, chem)        
        retcode = gsa_problem(params, prob, t_span, gsa_p, soln_cntrl)               
        return retcode
    else
        println("Models yet to be implemented\n")
    end
end
    

"""
A common function for all reactor models to define the sensitivity problem 
#   Usage 
    gsa_problem(params, prob, gsa_p)
-   params : a tuple of parameters returned by the respective reactor model invoked 
-   prob : the problem definition for the reactor model 
-   gsa_p : struct of the type gsa_parameters
"""
function gsa_problem(params, prob, t_span, gsa_p::Gsa_parameters, soln_cntrl::Solver)
    gsa_smech_params = Array{Float64,1}()
    state, thermo_obj, md, cp, chem = params        
    n_surface_species = length(md.sm.species) # number of surface species 
    n_gaspahse_species = length(thermo_obj.thermo_all) # number of gasphase species 
    n_species = n_gaspahse_species + n_surface_species # total number of species
    # get the parameters from the surface reaction mechaism  
    gsa_smech_parameter_map!(gsa_smech_params, gsa_p.gsa_srxn_ids, gsa_p.gsa_srxn_constraint_ids, md)
    println("Sample size: ", gsa_p.N)
    println("Number of parameters: ", length(gsa_p.gsa_srxn_ids))
    # create the lower limit for the parameters 
    lb = gsa_smech_params .* (1-0.01*gsa_p.lower_pc)
    # create the upper limit for the parameters 
    ub = gsa_smech_params .* (1+0.01*gsa_p.upper_pc)    
    
    parameter_ratio = Array{Float64,1}()  # Array for storing the original parameter ratio (forward to reverse)
    if count(x->x>0, gsa_p.gsa_srxn_constraint_ids) > 0
        # get the initial parameter ratios if the reverse reaction constraints are specified
        initial_parameter_ratio!(parameter_ratio, gsa_p.gsa_srxn_ids, gsa_p.gsa_srxn_constraint_ids,md)          
    end 
    
    retcode = perform_gsa(params,prob, t_span, gsa_p, lb, ub, parameter_ratio, soln_cntrl)
    return retcode

end


function perform_gsa(params, prob, t_span, gsa_p::Gsa_parameters,lb::Array{Float64,1}, ub::Array{Float64,1}, parameter_ratio::Array{Float64,1}, soln_cntrl::Solver)
    state, thermo_obj, md, cp, chem = params 

    #Get the name of all gasphase species 
    species_list = Array{String,1}()
    for sp_thermo in thermo_obj.thermo_all
        push!(species_list,sp_thermo.name)
    end
    append!(species_list,md.sm.species)

    # println(species_list)

    #get the id of species to be monitored from the gsa_species list
    monitor_species_ids = map(x->get_index(x,species_list),gsa_p.monitor)    

    
    
    func_eval_count = 0
    t_calls = gsa_p.N*(2+length(gsa_p.gsa_srxn_ids))         
    println("Total function calls :", t_calls)
    cb = FunctionCallingCallback(monitor_integration)
    model_response = zeros(length(monitor_species_ids))
    
    function sens(gsa_params::Array{Float64})                        
        func_eval_count += 1
        println(func_eval_count-t_calls)
        # println("Parameters ", gsa_params)
        SurfaceReactions.update_params!(md,gsa_params,gsa_p.gsa_srxn_ids,gsa_p.gsa_srxn_constraint_ids,parameter_ratio)
        params_updated = (state, thermo_obj, md, cp, chem)
        s_prob = remake(prob,p=params_updated)                
        if lowercase(strip(gsa_model)) == "plug"
            sol = solve(s_prob, CVODE_BDF(), reltol=soln_cntrl.reltol, abstol=soln_cntrl.abstol, save_everystep=false, callback=cb)
            for k in eachindex(monitor_species_ids)                
                model_response[k] = sol.u[end][monitor_species_ids[k]]
            end
        else
            sol = solve(SteadyStateProblem(s_prob), DynamicSS(CVODE_BDF()), dt=1e-10, reltol=soln_cntrl.reltol, abstol=soln_cntrl.abstol, save_everystep=false)
            for k in eachindex(monitor_species_ids)                
                model_response[k] = sol.u[monitor_species_ids[k]]
            end
        end
        
        return model_response
    end


    sampler = SobolSample()    
    A,B = QuasiMonteCarlo.generate_design_matrices(gsa_p.N,lb,ub,sampler)

    sobol_indices = gsa(sens,Sobol(order=[0,1]),A,B)

    
    sobol_first, sobol_total = o_streams  
    
    
    for k in 1:length(gsa_p.gsa_srxn_ids)
        data_out1 = Array{Float64,1}()
        data_out2 = Array{Float64,1}()
        for i in eachindex(monitor_species_ids)
            push!(data_out1, sobol_indices.S1[i,k])
            push!(data_out2, sobol_indices.ST[i,k])
        end
        p = "k"*string(gsa_p.gsa_srxn_ids[k])        
        write_out(sobol_first,  p,  data_out1)
        write_out(sobol_total,  p,  data_out2)
    end
    close(sobol_first)
    close(sobol_total)
    return :Success
    # # output file
    # sobol_total = open("sobol_total.dat","w")     
    # sobol_first = open("sobol_first.dat","w")
    # sobol_second = open("sobol_second.dat","w")
    # header_data = [ "p"*string(k) for k in gsa_p.gsa_srxn_ids]
    # create_header(sobol_total,["species"],header_data)
    # create_header(sobol_first,["species"],header_data)
    # create_header(sobol_second,["species"],header_data)

    # for i in eachindex(monitor_species_ids)
    #     println(species_list[i])
    #     write_out(sobol_total, species_list[i], gsa_p.monitor[i],sobol_indices.ST[i,:])
    #     write_out(sobol_first, species_list[i], gsa_p.monitor[i],sobol_indices.ST[i,:])
    #     write_out(sobol_second, species_list[i], gsa_p.monitor[i],sobol_indices.ST[i,:])
    # end    

    # for i in eachindex(header_data)       
    #     write_out(sobol_total, header_data[i], gsa_p.monitor[i],sobol_indices.ST[i,:])
    #     write_out(sobol_first, species_list[i], gsa_p.monitor[i],sobol_indices.ST[i,:])
        # write_out(sobol_second, species_list[i], gsa_p.monitor[i],sobol_indices.ST[i,:])
    # end

    # close(sobol_total)
    # close(sobol_first)
    # close(sobol_second)


    # function p_matrix!(m_resp, PM)        
    #     n_rows = size(PM,1)                
    #     #=
    #     each row will represent a full set of parameters for the model.
    #     Calculate the model response for all parameter space of PM. For 
    #     every set of parameters the code will return response of the 
    #     model for all species that are being monitored
    #     =#
    #     for i in 1:n_rows
    #         gsap = PM[i,:]                     
    #         # mres = sens(gsap)              
    #         mres = ishigami(gsap)
    #         # println(gsap, "\t", mres)   
    #         push!(m_resp, mres)                  
    #     end        
    # end
    
    # println("Evaluating model for A matrix parameters")
    # yA = Array{Array{Float64,1},1}()
    # p_matrix!(yA, A)    
    # display(yA)
    # println("Evaluating model for B matrix parameters")
    # yB = Array{Array{Float64,1},1}()
    # p_matrix!(yB, B)
    # display(yB)
    # #=
    # replace i'th column of A with the i'th column of B.
    # This replacement procedure will create k matrices, where k 
    # is the number of parameters 
    # =#
    # println("Evaluating model for AB matrix parameters")
    # yAB = Dict{Int64, Array{Array{Float64,1},1}}()
    # res_AB = Array{Array{Float64,1},1}()
    # n_col = size(A,2)
    # for i in 1:n_col
    #     C = A[:,i]
    #     A[:,i] = B[:,i] 
    #     println("_______________")
    #     # display(A)       
    #     p_matrix!(res_AB, A)
    #     yAB[i] = deepcopy(res_AB)
    #     display(yAB[i])
    #     empty!(res_AB)
    #     A[:,i] = C       
    # end
    

    # println("Model response evaluation for parameter space completed")
    # sp_id = 1
    # S1 = Array{Float64,1}()
    # ST = Array{Float64,1}()
    # for i in 1:3
    #     x = y = z = q = 0.0
    #     for j in 1:gsa_p.N
    #         x += yB[j][sp_id]*( yAB[i][j][sp_id] - yA[j][sp_id] )
    #         y += (yA[j][sp_id])^2
    #         z += yA[j][sp_id]                
    #         q += (yA[j][sp_id] - yAB[i][j][sp_id])^2
    #     end
    #     Dr = (y/gsa_p.N) - (z/gsa_p.N)^2
    #     push!(S1, (x/gsa_p.N)/Dr)
    #     push!(ST, q/(2*gsa_p.N)/Dr)        
    # end
    # println(S1)
    # println(ST)
    
    # header_data = [ "p"*string(k) for k in gsa_p.gsa_srxn_ids]
    # sp_id = 0
    # for m in gsa_p.monitor        
    #     # fname = "gsa_first_"*m*".dat"
    #     # tname = "gsa_total_"*m*".dat"
    #     # f_index = open(fname, "w")
    #     # t_index = open(tname, "w")
    #     # create_header(f_index, header_data)
    #     # create_header(t_index, header_data)
    #     S1 = Array{Float64,1}()
    #     ST = Array{Float64,1}()
    #     sp_id += 1
    #     for i in 1:length(gsa_p.gsa_srxn_ids)
    #         x = y = z = q = 0.0                        
    #         for j in 1:gsa_p.N 
    #             x += yB[j][sp_id]*( yAB[i][j][sp_id] - yA[j][sp_id] )
    #             y += (yA[j][sp_id])^2
    #             z += yA[j][sp_id]                
    #             q += (yA[j][sp_id] - yAB[i][j][sp_id])^2
    #         end
    #         Dr = y/gsa_p.N - (z/gsa_p.N)^2
    #         push!(S1, (x/gsa_p.N)/Dr)
    #         push!(ST, q/(2*gsa_p.N)/Dr)
    #     end
    #     println(S1)
    #     println(ST)
    #     # write_out(f_index, S1)
    #     # write_out(t_index, ST)
    #     # empty!(S1)
    #     # empty!(ST)
    #     # close(f_index)
    #     # close(t_index)
    # end
    
    
end

function write_out(file_stream, args...)
    
    for i in eachindex(args)
        if isa(args[i], String)            
            @printf(file_stream,"%10s\t",args[i])
        else            
            for k in args[i]
                @printf(file_stream,"%+.4e\t",k)
            end
        end
    end
    @printf(file_stream,"\n")
end

function monitor_integration(u, t, integrator)
    @printf("%.4e\n",t)
end
"""
A function to extract the parameters for sensitivity analysis from the xml input 
    The function returns two integer arrays. forward_rxn_ids contains the ids of the 
    reactions whoes parameters will be considered as sensitivity parameters.
    reverse_rxn_ids contains the ids of the reactions aginst with the forward reaction 
    parameters are correlated. 
"""
function get_rxn_ids(xmlroot::XMLElement, xmltag::AbstractString)
    N = 10000
    p_mech = get_elements_by_tagname(xmlroot,xmltag)
    p_mech_content = ""
    if is_elementnode(p_mech[1])
        p_mech_content = content(p_mech[1])        
        
        p_mech_content_split = split(p_mech_content,",")           
        forward_rxn_ids = Array{Int64,1}()
        reverse_rxn_ids = Array{Int64,1}()
        if strip(p_mech_content_split[1]) == ":"        
            return [i for i in 1:N], [0 for i in 1:N]
        end
        for items in p_mech_content_split
            if occursin("=",items)
                rxn_sets = split(items,"=")        
                forward_sets = rxn_sets[1]
                reverse_sets = rxn_sets[2]
                forward_sets = map(x->parse(Int64,x),split(forward_sets,":"))
                reverse_sets = map(x->parse(Int64,x),split(reverse_sets,":"))
                if length(forward_sets) > 1
                    for ids in forward_sets[1]:forward_sets[2]
                        push!(forward_rxn_ids,ids)
                    end
                    for ids in reverse_sets[1]:reverse_sets[2]
                        push!(reverse_rxn_ids,ids)
                    end
                else                
                    push!(forward_rxn_ids,forward_sets[1])
                    push!(reverse_rxn_ids,reverse_sets[1])
                end       
            else            
                forward_sets = map(x->parse(Int64,x),split(items,":"))            
                if length(forward_sets) > 1
                    for ids in forward_sets[1]:forward_sets[2]
                        push!(forward_rxn_ids,ids)
                        push!(reverse_rxn_ids,0)
                    end
                else
                    push!(forward_rxn_ids,forward_sets[1])
                    push!(reverse_rxn_ids,0)
                end
            end
        end
        return (forward_rxn_ids, reverse_rxn_ids)
    else
        nothing, nothing
    end
end

"""
This function returns the sticking coefficient or the pre-exponent factor 
    corresponding to the reaction ids present in gsa_rxn_ids 
#   Usage
    gsa_smech_parameter_map!(gsa_smech_params, gsa_srxn_ids, gsa_srxn_constraint_ids, md )    
-   gsa_smech_params : vector containing the parameters extracted from the inoput mechanism (output)
-   gsa_srxn_ids : reaction ids whoes parameters needs to be extracted from the mechanism (input)
-   gsa_srxn_constraint_ids : id of the reactions which needs to be considered for reverse reaction 
-   md: MechanismDefinition
"""
function gsa_smech_parameter_map!(gsa_smech_params::Array{Float64,1}, gsa_srxn_ids::Array{Int64,1}, gsa_srxn_constraint_ids::Array{Int64,1}, md )    
    
    function populate(rxn)
        if rxn.type == ReactionCommons.stick
            push!(gsa_smech_params,rxn.params.s)
        else
            push!(gsa_smech_params,rxn.params.k0)
        end
    end
    # if the length of gsa_srxn_ids is very high then it means all reactions parameters need to be considered
    if length(gsa_srxn_ids) >= length(md.sm.reactions)
        for rxn in md.sm.reactions
            populate(rxn)
        end
        gsa_srxn_ids_copy = copy(gsa_srxn_ids)                
        empty!(gsa_srxn_ids)
        gsa_srxn_constraint_ids_copy = copy(gsa_srxn_constraint_ids)        
        empty!(gsa_srxn_constraint_ids)
        for i in 1:length(md.sm.reactions)
            push!(gsa_srxn_ids,gsa_srxn_ids_copy[i])
            push!(gsa_srxn_constraint_ids,gsa_srxn_constraint_ids_copy[i])
        end
        empty!(gsa_srxn_ids_copy)
        empty!(gsa_srxn_constraint_ids_copy)
    else
        for rxn in md.sm.reactions
            if rxn.id in gsa_srxn_ids
                populate(rxn)
            end        
        end
    end

    return gsa_smech_params
end


"""
A function to calculate the intial ratio of pre-exponential factors
"""
function initial_parameter_ratio!(p_ratio::Array{Float64,1},rxn_ids::Array{Int64,1},constraints_ids::Array{Int64,1}, md)
    for i in eachindex(rxn_ids)        
        if constraints_ids[i] != 0
            if md.sm.reactions[rxn_ids[i]].type == ReactionCommons.stick
                ratio = md.sm.reactions[rxn_ids[i]].params.k0/md.sm.reactions[constraints_ids[i]].params.k0                
                push!(p_ratio,ratio)
            else                
                ratio = md.sm.reactions[rxn_ids[i]].params.k0/md.sm.reactions[constraints_ids[i]].params.k0
                push!(p_ratio,ratio)
            end
        else
            push!(p_ratio,0.0)
        end
    end
end

#end of module
end
