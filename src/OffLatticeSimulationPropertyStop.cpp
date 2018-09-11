#include "../../../projects/ISP/src/OffLatticeSimulationPropertyStop.hpp"

#include <boost/make_shared.hpp>

#include "CellBasedEventHandler.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "StepSizeException.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
OffLatticeSimulationPropertyStop<ELEMENT_DIM,SPACE_DIM>::OffLatticeSimulationPropertyStop(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                                bool deleteCellPopulationInDestructor,
                                                bool initialiseCells
                                                )
    : AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>(rCellPopulation, deleteCellPopulationInDestructor, initialiseCells),
    p_property()
{
    if (!dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation))
    {
        EXCEPTION("OffLatticeSimulations require a subclass of AbstractOffLatticeCellPopulation.");
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool OffLatticeSimulationPropertyStop<ELEMENT_DIM,SPACE_DIM>::StoppingEventHasOccurred()
{
    if(p_property->GetCellCount()<1){
		return true;
	}
	else{
		return false;
	}
}


//Public access to Stopping Event bool
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool OffLatticeSimulationPropertyStop<ELEMENT_DIM,SPACE_DIM>::HasStoppingEventOccurred()
{
    return StoppingEventHasOccurred();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulationPropertyStop<ELEMENT_DIM,SPACE_DIM>::SetStopProperty(boost::shared_ptr<AbstractCellProperty> stopPropertySetting)
{
	 p_property = stopPropertySetting;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulationPropertyStop<ELEMENT_DIM,SPACE_DIM>::AddForce(boost::shared_ptr<AbstractForce<ELEMENT_DIM,SPACE_DIM> > pForce)
{
    mForceCollection.push_back(pForce);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulationPropertyStop<ELEMENT_DIM,SPACE_DIM>::RemoveAllForces()
{
    mForceCollection.clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulationPropertyStop<ELEMENT_DIM,SPACE_DIM>::AddCellPopulationBoundaryCondition(boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > pBoundaryCondition)
{
    mBoundaryConditions.push_back(pBoundaryCondition);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulationPropertyStop<ELEMENT_DIM,SPACE_DIM>::RemoveAllCellPopulationBoundaryConditions()
{
    mBoundaryConditions.clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulationPropertyStop<ELEMENT_DIM,SPACE_DIM>::SetNumericalMethod(boost::shared_ptr<AbstractNumericalMethod<ELEMENT_DIM, SPACE_DIM> > pNumericalMethod)
{
    mpNumericalMethod = pNumericalMethod;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const boost::shared_ptr<AbstractNumericalMethod<ELEMENT_DIM, SPACE_DIM> > OffLatticeSimulationPropertyStop<ELEMENT_DIM,SPACE_DIM>::GetNumericalMethod() const
{
    return mpNumericalMethod;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >& OffLatticeSimulationPropertyStop<ELEMENT_DIM,SPACE_DIM>::rGetForceCollection() const
{
    return mForceCollection;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulationPropertyStop<ELEMENT_DIM,SPACE_DIM>::UpdateCellLocationsAndTopology()
{
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::POSITION);

    double time_advanced_so_far = 0;
    double target_time_step  = this->mDt;
    double present_time_step = this->mDt;

    while (time_advanced_so_far < target_time_step)
    {
        // Store the initial node positions (these may be needed when applying boundary conditions)
        std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > old_node_locations;

        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
             node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            old_node_locations[&(*node_iter)] = (node_iter)->rGetLocation();
        }

        // Try to update node positions according to the numerical method
        try
        {
            mpNumericalMethod->UpdateAllNodePositions(present_time_step);
            ApplyBoundaries(old_node_locations);

            // Successful time step! Update time_advanced_so_far
            time_advanced_so_far += present_time_step;

            // If using adaptive timestep, then increase the present_time_step (by 1% for now)
            if (mpNumericalMethod->HasAdaptiveTimestep())
            {
                ///\todo #2087 Make this a settable member variable
                double timestep_increase = 0.01;
                present_time_step = std::min((1+timestep_increase)*present_time_step, target_time_step - time_advanced_so_far);
            }

        }
        catch (StepSizeException& e)
        {
            // Detects if a node has travelled too far in a single time step
            if (mpNumericalMethod->HasAdaptiveTimestep())
            {
                // If adaptivity is switched on, revert node locations and choose a suitably smaller time step
                RevertToOldLocations(old_node_locations);
                present_time_step = std::min(e.GetSuggestedNewStep(), target_time_step - time_advanced_so_far);
            }
            else
            {
                // If adaptivity is switched off, terminate with an error
                EXCEPTION(e.what());
            }
        }
    }

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::POSITION);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulationPropertyStop<ELEMENT_DIM,SPACE_DIM>::RevertToOldLocations(std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > oldNodeLoctions)
{
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
        ++node_iter)
    {
        (node_iter)->rGetModifiableLocation() = oldNodeLoctions[&(*node_iter)];
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulationPropertyStop<ELEMENT_DIM,SPACE_DIM>::ApplyBoundaries(std::map<Node<SPACE_DIM>*,c_vector<double, SPACE_DIM> > oldNodeLoctions)
{
    // Apply any boundary conditions
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > >::iterator bcs_iter = mBoundaryConditions.begin();
         bcs_iter != mBoundaryConditions.end();
         ++bcs_iter)
    {
        (*bcs_iter)->ImposeBoundaryCondition(oldNodeLoctions);
    }

    // Verify that each boundary condition is now satisfied
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > >::iterator bcs_iter = mBoundaryConditions.begin();
         bcs_iter != mBoundaryConditions.end();
         ++bcs_iter)
    {
        if (!((*bcs_iter)->VerifyBoundaryCondition()))
        {
            EXCEPTION("The cell population boundary conditions are incompatible.");
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulationPropertyStop<ELEMENT_DIM,SPACE_DIM>::WriteVisualizerSetupFile()
{
    if (PetscTools::AmMaster())
    {
        for (unsigned i=0; i<this->mForceCollection.size(); i++)
        {
            this->mForceCollection[i]->WriteDataToVisualizerSetupFile(this->mpVizSetupFile);
        }

        this->mrCellPopulation.WriteDataToVisualizerSetupFile(this->mpVizSetupFile);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulationPropertyStop<ELEMENT_DIM,SPACE_DIM>::SetupSolve()
{
    // Clear all forces
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        node_iter->ClearAppliedForce();
    }

    // Use a forward Euler method by default, unless a numerical method has been specified already
    if (mpNumericalMethod == nullptr)
    {
        mpNumericalMethod = boost::make_shared<ForwardEulerNumericalMethod<ELEMENT_DIM, SPACE_DIM> >();
    }
    mpNumericalMethod->SetCellPopulation(dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(this->mrCellPopulation)));
    mpNumericalMethod->SetForceCollection(&mForceCollection);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulationPropertyStop<ELEMENT_DIM,SPACE_DIM>::OutputAdditionalSimulationSetup(out_stream& rParamsFile)
{
    // Loop over forces
    *rParamsFile << "\n\t<Forces>\n";
    for (typename std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM,SPACE_DIM> > >::iterator iter = mForceCollection.begin();
         iter != mForceCollection.end();
         ++iter)
    {
        // Output force details
        (*iter)->OutputForceInfo(rParamsFile);
    }
    *rParamsFile << "\t</Forces>\n";

    // Loop over cell population boundary conditions
    *rParamsFile << "\n\t<CellPopulationBoundaryConditions>\n";
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > >::iterator iter = mBoundaryConditions.begin();
         iter != mBoundaryConditions.end();
         ++iter)
    {
        // Output cell boundary condition details
        (*iter)->OutputCellPopulationBoundaryConditionInfo(rParamsFile);
    }
    *rParamsFile << "\t</CellPopulationBoundaryConditions>\n";

    // Output numerical method details
    *rParamsFile << "\n\t<NumericalMethod>\n";
    mpNumericalMethod->OutputNumericalMethodInfo(rParamsFile);
    *rParamsFile << "\t</NumericalMethod>\n";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulationPropertyStop<ELEMENT_DIM,SPACE_DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::OutputSimulationParameters(rParamsFile);
}

// Explicit instantiation
template class OffLatticeSimulationPropertyStop<1,1>;
template class OffLatticeSimulationPropertyStop<1,2>;
template class OffLatticeSimulationPropertyStop<2,2>;
template class OffLatticeSimulationPropertyStop<1,3>;
template class OffLatticeSimulationPropertyStop<2,3>;
template class OffLatticeSimulationPropertyStop<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(OffLatticeSimulationPropertyStop)
