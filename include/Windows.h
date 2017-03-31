
#include <cstddef> //For size_t
#include <vector>

namespace myFuncs
{  

	//-------------------------------------------------------------
	//Base window class
	//-------------------------------------------------------------
	class WindowBase {
		
	public:
		
		//Constructor
		WindowBase(const size_t windowLength);
				
		//Eval window at index
		virtual double eval(const long index) const = 0;
		
		//Get window length
		inline size_t getWindowLength() const {return m_windowLength;}
		
		//Get window values
		//NOTE - the size of the vector is getWindowLength() + 1, so that the last value is 0
		std::vector<double> getWindowValues() const;
		
				
		template <typename T>
		std::vector<double> windowAnInput(const std::vector<T>& inputs) const {
			std::vector<double> outputs;
			outputs.reserve(inputs.size());
			
			for(size_t index = 0; index < inputs.size(); ++index)
				outputs.push_back(inputs[index] * eval(index) );
			
			return outputs;
		}
		
		
	private:
		//Size of window
		size_t m_windowLength;
		
	};
	
	class Hamming : public WindowBase {
		
	public:
		Hamming(const size_t windowLength);
		
		double eval(const long index) const override final;
	};
	
	class Hann : public WindowBase {
	public:
		Hann(const size_t windowLength);
		
		double eval(const long index) const override final;
		
	};	
	
}